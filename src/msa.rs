use std::{
    cmp,
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    path::Path,
};

use noodles::sam::alignment::RecordBuf;

use crate::{
    cli::Cli,
    kits::{
        bam_ext::{BamRecordExt, aligned_pairs_full},
        pbar::{DEFAULT_INTERVAL, get_bar_pb, get_spin_pb},
    },
};

pub fn msa_entrance(cli: &Cli) {
    let filepath = cli.metric_args.io_args.query[0].clone();
    // 等价于原 rust-htslib `set_threads(cli.threads.unwrap_or(num_cpus))`:
    // BGZF 解压 worker 数。
    let file = fs::File::open(&filepath).unwrap();
    let worker_count =
        std::num::NonZero::new(cli.threads.unwrap_or(num_cpus::get()).max(1)).unwrap();
    let bgzf_reader = noodles::bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = noodles::bam::io::Reader::from(bgzf_reader);
    let header = reader.read_header().unwrap();
    let output_filepath = cli.metric_args.get_oup_file(cli.mode);
    // 用 Path 换扩展名（bug-report.md Bug 9：旧 rsplit_once 依赖路径含「.」，
    // 无扩展名时 panic，目录名含「.」时会切错位置）
    let output_filepath = Path::new(&output_filepath)
        .with_extension("fasta")
        .to_string_lossy()
        .into_owned();

    let shannon_entroy_output_filepath = format!("{}.shannon.txt", output_filepath);

    let pb = get_spin_pb(format!("reading {}", filepath), DEFAULT_INTERVAL);
    let all_records = reader
        .record_bufs(&header)
        .inspect(|_rec| {
            pb.inc(1);
        })
        .flatten()
        .filter(|rec| !rec.flags().is_unmapped())
        .filter(|rec| !rec.flags().is_secondary())
        .filter(|rec| !rec.flags().is_supplementary())
        .collect::<Vec<_>>();
    pb.finish();

    let major_pos_ins = compute_max_ins_of_each_ref_position(&all_records, Some(0), None);

    let mut major_pos_ins_vec = major_pos_ins
        .iter()
        .map(|(&k, &v)| (k as usize, v as usize))
        .collect::<Vec<_>>();
    major_pos_ins_vec.sort_by_key(|v| v.0);

    let major = major_pos_ins_vec
        .iter()
        .flat_map(|&(major_pos, ins_size)| vec![major_pos; ins_size + 1].into_iter())
        .collect::<Vec<_>>();

    let mut cursor = 0;
    let major_start_point = major_pos_ins_vec
        .iter()
        .map(|&(_, max_ins)| {
            let cur_point = cursor;
            cursor += max_ins + 1;
            cur_point
        })
        .collect::<Vec<_>>();
    let major_pos2major_starting_point = major_pos_ins_vec
        .iter()
        .map(|&(ma, _)| ma)
        .zip(major_start_point)
        .collect::<HashMap<_, _>>();

    let msa_len = major.len();
    draw_msa(
        &all_records,
        0,
        usize::MAX,
        &major_pos2major_starting_point,
        msa_len,
        &output_filepath,
        &shannon_entroy_output_filepath,
    );
}

fn compute_max_ins_of_each_ref_position(
    records: &Vec<RecordBuf>,
    rstart: Option<usize>,
    rend: Option<usize>,
) -> HashMap<i64, i32> {
    let mut pos2ins = HashMap::new();

    let rstart = rstart.map(|v| v as i64);
    let rend = rend.map(|v| v as i64);

    for record in records {
        let record_ext = BamRecordExt::new(record);
        let mut start = rstart.unwrap_or(record_ext.reference_start() as i64);
        let mut end = rend.unwrap_or(record_ext.reference_end() as i64);
        // start 侧取 max（下限夹到 record 起点），end 侧取 min（上限夹到
        // record 末端）。旧代码 end 侧误用 max，调用方传入的 rend 上限会被
        // 静默放大回 record 末尾（bug-report.md Bug 6；当前唯一调用点
        // rend=None 未触发）。
        start = start.max(record_ext.reference_start() as i64);
        end = end.min(record_ext.reference_end() as i64);

        let mut rpos_cursor = None;
        let mut cur_ins = 0;
        let query_end = record_ext.query_alignment_end();
        let mut aligned_pair_full = aligned_pairs_full(record).collect::<Vec<_>>();
        if aligned_pair_full.is_empty() {
            continue;
        }

        // this make the following for loop correct.
        // if the query match to the last base of the ref sequence. the following for loop won't give the right result
        // but add this , it will get the right result.
        if let Some(last_ref_pos) = aligned_pair_full.last().unwrap()[0] {
            aligned_pair_full.push([Some(last_ref_pos + 1), None]);
        }

        for [qpos, rpos] in aligned_pair_full.into_iter() {
            if rpos.is_some() {
                rpos_cursor = rpos;
            }
            if rpos_cursor.is_none() {
                continue;
            }

            if rpos_cursor.unwrap() < start {
                continue;
            }

            if rpos_cursor.unwrap() >= end {
                // set the last rpos max ins and then break
                let rpos_ = rpos_cursor.unwrap();
                pos2ins.entry(rpos_ - 1).or_insert(0);
                *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                    cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                break;
            }

            if let Some(qpos_) = qpos
                && qpos_ as usize >= query_end
            {
                // set the last rpos max ins and then break
                let rpos_ = rpos_cursor.unwrap();

                pos2ins.entry(rpos_).or_insert(0);
                *pos2ins.get_mut(&rpos_).unwrap() =
                    cmp::max(*pos2ins.get(&rpos_).unwrap(), cur_ins);
                break;
            }

            if let Some(rpos_) = rpos {
                if rpos_ > start {
                    pos2ins.entry(rpos_ - 1).or_insert(0);
                    *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                        cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                }

                cur_ins = 0;
            } else {
                cur_ins += 1;
            }
        }
    }

    pos2ins
}

#[derive(Debug, Clone)]
struct LocusBaseCounter {
    counter: HashMap<char, usize>,
}
impl Default for LocusBaseCounter {
    fn default() -> Self {
        let mut counter = HashMap::new();
        counter.insert('-', 0);
        counter.insert('A', 0);
        counter.insert('C', 0);
        counter.insert('G', 0);
        counter.insert('T', 0);
        Self { counter }
    }
}

impl LocusBaseCounter {
    fn shannon_entropy(&self) -> f32 {
        let mut total_count = 0;
        for &count in self.counter.values() {
            total_count += count;
        }

        if total_count == 0 {
            return 0.0;
        }

        let mut entropy = 0.0;
        for &count in self.counter.values() {
            if count > 0 {
                let p = count as f32 / total_count as f32;
                entropy -= p * p.log2();
            }
        }
        entropy
    }
}

fn draw_msa(
    records: &[RecordBuf],
    rstart: usize,
    rend: usize,
    major_pos2major_starting_point: &HashMap<usize, usize>,
    msa_len: usize,
    ofilename: &str,
    shann_entropy_ofilename: &str,
) {
    let mut all_locus_base_counter = vec![LocusBaseCounter::default(); msa_len];

    let ofile = fs::File::create(ofilename).unwrap();
    let mut writer = BufWriter::new(ofile);

    let mut single_record_msa = vec!['-'; msa_len];

    let pb = get_bar_pb(
        format!("dumping msa result to {}", ofilename),
        DEFAULT_INTERVAL,
        records.len() as u64,
    );

    for record in records.iter() {
        single_record_msa.fill('-');

        let record_ext = BamRecordExt::new(record);
        let rstart = record_ext.reference_start().max(rstart);
        let rend = record_ext.reference_end().min(rend);
        let qstart = record_ext.query_alignment_start();
        let qend = record_ext.query_alignment_end();

        let mut qpos_cursor = None;
        let mut rpos_cursor = None;

        let query_seq = record_ext.get_seq();

        let query_seq_bytes = query_seq.chars().collect::<Vec<_>>();
        let mut minor_cursor = 0;

        for [qpos, rpos] in aligned_pairs_full(record) {
            if rpos.is_some() {
                rpos_cursor = rpos;
            }
            if qpos.is_some() {
                qpos_cursor = qpos;
            }

            if rpos_cursor.is_none()
                || qpos_cursor.is_none()
                || (rpos_cursor.unwrap_or(0) as usize) < rstart
                || (qpos_cursor.unwrap_or(0) as usize) < qstart
            {
                continue;
            }

            if rpos_cursor.unwrap() as usize >= rend || qpos_cursor.unwrap() as usize >= qend {
                break;
            }

            if rpos.is_some() {
                minor_cursor = 0;
            } else {
                minor_cursor += 1;
            }

            // fill base and dw
            // get 不到该 ref 位点（覆盖它的 read 全部以 deletion 收尾且该位点
            // 是各自的 reference_end，从未进入 pos2ins）时跳过该柱，不再 panic
            // （bug-report.md Bug 6）。
            let Some(major_start) =
                major_pos2major_starting_point.get(&(rpos_cursor.unwrap() as usize))
            else {
                continue;
            };
            let tt = major_start + minor_cursor;
            if let Some(qpos_) = qpos {
                single_record_msa[tt] = query_seq_bytes[qpos_ as usize];
            }
        }

        single_record_msa
            .iter()
            .zip(all_locus_base_counter.iter_mut())
            .for_each(|(&base, base_counter)| {
                // or_insert(0)：非 -/A/C/G/T 碱基（如 N）也要计入分母
                // （bug-report.md Bug 10：旧 and_modify 只更新预置 key，N 被静默丢弃）
                *base_counter.counter.entry(base).or_insert(0) += 1;
            });

        writeln!(&mut writer, ">{}", record_ext.get_qname()).unwrap();
        writeln!(
            &mut writer,
            "{}",
            single_record_msa.iter().collect::<String>()
        )
        .unwrap();

        pb.inc(1);
    }
    pb.finish();

    // 显式 flush 并检查错误：drop-flush 会静默吞掉磁盘满等错误，导致
    // 「成功」退出却丢尾部数据
    writer
        .flush()
        .unwrap_or_else(|e| panic!("flush msa output {ofilename} failed: {e}"));

    let entropy_str = all_locus_base_counter
        .iter()
        .map(|locus| locus.shannon_entropy().to_string())
        .collect::<Vec<_>>()
        .join(",");

    let mut ofile = fs::File::create(shann_entropy_ofilename).unwrap();
    ofile.write_all(entropy_str.as_bytes()).unwrap();
    ofile.flush().unwrap();
}
