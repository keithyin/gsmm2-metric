#!/usr/bin/env python3
"""生成 bam2-msa / metric BAM 路径的小型 BAM fixture。

本机没有 samtools/pysam,这里按 BAM 规范 + BGZF(RFC1952 + BC extra field)
手工编码。reads 取自 test_data 参考序列的子串,并覆盖:
  EQX cigar(=/X)、soft clip、insertion、deletion、hard clip、reverse strand、
  secondary/supplementary/unmapped、qual 缺失(*)、np/rq 缺失,
  以及 np(int32/uint8 两种编码)、rq(float)、cx/ch 标量、be/dw(B 数组)tag。

用法: python3 scripts/gen_bam_fixture.py [out.bam]
"""
import re
import struct
import sys
import zlib

REF = "test_data/ref_Saureus_ATCC25923.m.new.corrected.fasta"

CIGAR_OPS = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
# BAM 规范的 4-bit 编码(IUPAC16): =0 A=1 C=2 M=3 G=4 ... T=8 ... N=15
SEQ_CODES = {"A": 1, "C": 2, "G": 4, "T": 8, "N": 15}


def reg2bin(beg, end):
    end -= 1
    if beg >> 14 == end >> 14:
        return 4681 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return 585 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return 73 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return 9 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return 1 + (beg >> 26)
    return 0


def pack_seq(seq):
    seq = seq.upper()
    out = bytearray()
    for i in range(0, len(seq), 2):
        hi = SEQ_CODES[seq[i]]
        lo = SEQ_CODES[seq[i + 1]] if i + 1 < len(seq) else 0
        out.append((hi << 4) | lo)
    return bytes(out)


def parse_cigar(cigar):
    return [(int(n), CIGAR_OPS[op]) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def ref_span(cigar):
    return sum(n for n, op in parse_cigar(cigar) if op in (0, 2, 3, 7, 8))


def aux_i32(tag, v):
    return tag + b"i" + struct.pack("<i", v)


def aux_u8(tag, v):
    return tag + b"C" + struct.pack("<B", v)


def aux_u16(tag, v):
    return tag + b"S" + struct.pack("<H", v)


def aux_f32(tag, v):
    return tag + b"f" + struct.pack("<f", v)


def aux_barr(tag, subtype, values):
    fmt = {"C": "<B", "I": "<I", "S": "<H", "f": "<f"}[subtype]
    body = b"".join(struct.pack(fmt, v) for v in values)
    return tag + b"B" + subtype.encode() + struct.pack("<i", len(values)) + body


def aux_z(tag, s):
    return tag + b"Z" + s.encode() + b"\x00"


def bam_record(name, flag, ref_id, pos, mapq, cigar, seq, qual, aux=b""):
    """qual: phred 字节序列,或 None 表示 '*'。"""
    n_cigar_op = len(parse_cigar(cigar))
    l_seq = len(seq)
    end = pos + ref_span(cigar) if pos >= 0 else pos + 1
    bin_ = reg2bin(pos, end) if pos >= 0 else 4692
    if qual is None:
        quals = b"\xff" * l_seq
    else:
        assert len(qual) == l_seq
        quals = bytes(qual)
    body = (
        struct.pack("<i", ref_id)
        + struct.pack("<i", pos)
        + struct.pack("<BBH", len(name) + 1, mapq, bin_)
        + struct.pack("<HH", n_cigar_op, flag)
        + struct.pack("<iiii", l_seq, -1, -1, 0)
        + name.encode() + b"\x00"
        + struct.pack(f"<{n_cigar_op}I", *[n << 4 | op for n, op in parse_cigar(cigar)])
        + pack_seq(seq)
        + quals
        + aux
    )
    return struct.pack("<i", len(body)) + body


def bgzf_block(payload):
    co = zlib.compressobj(6, zlib.DEFLATED, -15)
    data = co.compress(payload) + co.flush()
    bsize = 26 + len(data) - 1  # header(18) + data + crc(4) + isize(4), -1 per spec
    header = struct.pack(
        "<BBBBIBBHBBH",
        0x1F, 0x8B, 0x08, 0x04, 0, 0, 255, 6,
        ord("B"), ord("C"), 2,
    ) + struct.pack("<H", bsize)
    return header + data + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def read_ref():
    with open(REF) as f:
        header = f.readline().strip().lstrip(">")
        seq = "".join(line.strip() for line in f)
    return header, seq


def main(out_path):
    ref_name, ref = read_ref()
    assert len(ref) > 1000

    def sub(a, b):
        return ref[a:b]

    # (name, flag, pos, mapq, cigar, seq, qual_phred, aux)
    # 除标注外 qual = 30;seq 与 cigar 相对参考保持一致(bam2-msa 直接用 cigar)。
    # reads 取 500bp(与 queries.fa 的合成 reads 同长度,保证 bandmap 能比对出 hits)。
    reads = []

    # r1: 完全匹配,np=i32
    s = sub(1000, 1500)
    reads.append(("r1_exact", 0, 1000, 60, "500=", s, 30,
                  aux_i32(b"np", 1) + aux_f32(b"rq", 0.95) + aux_u8(b"cx", 0)
                  + aux_u16(b"ch", 10) + aux_barr(b"be", "I", [1, 2, 3])
                  + aux_barr(b"dw", "C", [4, 5])))

    # r2: soft clip + 2 mismatch,np=i32,rq 略低(过滤用)
    s = sub(2000, 2500)
    mm1 = "A" if s[200] != "A" else "T"
    mm2 = "G" if s[201] != "G" else "C"
    q = "ACGTA" + s[:200] + mm1 + mm2 + s[202:500] + "TTTTT"
    reads.append(("r2_softclip_mm", 0, 2000, 60, "5S200=2X298=5S", q, 30,
                  aux_i32(b"np", 2) + aux_f32(b"rq", 0.90) + aux_u8(b"cx", 1)))

    # r3: insertion,np 用 uint8 编码
    s = sub(3000, 3500)
    q = s[:250] + "GT" + s[250:]
    reads.append(("r3_ins", 0, 3000, 60, "250=2I250=", q, 30,
                  aux_u8(b"np", 2) + aux_f32(b"rq", 0.85) + aux_z(b"zz", "hello")))

    # r4: deletion
    s = sub(4000, 4503)
    q = s[:250] + s[253:]
    reads.append(("r4_del", 0, 4000, 60, "250=3D250=", q, 30,
                  aux_i32(b"np", 3) + aux_f32(b"rq", 0.95)))

    # r5: 末尾 trailing insertion(bam2-msa 边界补丁路径)
    s = sub(5000, 5500)
    q = s + "AA"
    reads.append(("r5_trailing_ins", 0, 5000, 60, "500=2I", q, 30,
                  aux_i32(b"np", 1) + aux_f32(b"rq", 0.70)))

    # r6: hard clip(aligned_pairs_full 应跳过)
    s = sub(6000, 6500)
    reads.append(("r6_hardclip", 0, 6000, 60, "5H500=5H", s, 30,
                  aux_i32(b"np", 2) + aux_f32(b"rq", 0.99)))

    # r7: reverse strand
    s = sub(7000, 7500)
    mm = "A" if s[250] != "A" else "T"
    q = s[:250] + mm + s[251:500] + "TTTTT"
    reads.append(("r7_reverse", 16, 7000, 60, "250=1X249=5S", q, 30,
                  aux_i32(b"np", 1) + aux_f32(b"rq", 0.95)))

    # r8: secondary(会被 metric 丢弃,bam2-msa 过滤)
    s = sub(8000, 8500)
    reads.append(("r8_secondary", 256, 8000, 3, "500=", s, 30,
                  aux_i32(b"np", 1) + aux_f32(b"rq", 0.95)))

    # r9: supplementary
    s = sub(9000, 9500)
    reads.append(("r9_supplementary", 2048, 9000, 60, "500=", s, 30,
                  aux_i32(b"np", 1) + aux_f32(b"rq", 0.95)))

    # r10: unmapped,无 np/rq
    reads.append(("r10_unmapped", 4, -1, 0, "500=", sub(10000, 10500), 30, b""))

    # r11: 无 qual(*),np 缺失
    s = sub(11000, 11500)
    reads.append(("r11_noqual", 0, 11000, 60, "500=", s, None,
                  aux_f32(b"rq", 0.95)))

    # r12: rq 明显低于 0.9(过滤用)
    s = sub(12000, 12500)
    reads.append(("r12_lowrq", 0, 12000, 60, "500=", s, 30,
                  aux_i32(b"np", 3) + aux_f32(b"rq", 0.75)))

    text = (
        "@HD\tVN:1.6\tSO:coordinate\n"
        f"@SQ\tSN:{ref_name}\tLN:{len(ref)}\n"
        "@PG\tID:gsmm2-fixture\tPN:gen_bam_fixture.py\n"
    ).encode() + b"\x00"

    payload = b"BAM\x01" + struct.pack("<i", len(text)) + text
    payload += struct.pack("<i", 1)
    payload += struct.pack("<i", len(ref_name) + 1) + ref_name.encode() + b"\x00"
    payload += struct.pack("<i", len(ref))

    records = []
    for name, flag, pos, mapq, cigar, seq, qual, aux in reads:
        ref_id = 0 if pos >= 0 else -1
        phred = None if qual is None else bytes([qual] * len(seq))
        records.append(bam_record(name, flag, ref_id, pos, mapq, cigar, seq, phred, aux))

    with open(out_path, "wb") as f:
        f.write(bgzf_block(payload + b"".join(records)))
        f.write(BGZF_EOF)

    print(f"wrote {out_path}: {len(reads)} records")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "test_data/queries_bam/fixture.bam")
