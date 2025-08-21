use std::{
    cmp,
    collections::HashMap,
    fs,
    io::{BufWriter, Write},
    usize,
};

use gskits::{
    gsbam::bam_record_ext::BamRecordExt,
    pbar::{DEFAULT_INTERVAL, get_bar_pb, get_spin_pb},
};
use rust_htslib::{
    self,
    bam::{Read, Record, ext::BamRecordExtensions},
};

use crate::cli::Cli;
pub fn msa_entrence(cli: &Cli) {
    let filepath = cli.metric_args.io_args.query[0].clone();
    let mut reader = rust_htslib::bam::Reader::from_path(&filepath).unwrap();
    reader
        .set_threads(cli.threads.unwrap_or(num_cpus::get()))
        .unwrap();
    let output_filepath = cli.metric_args.get_oup_file(cli.mode);
    let output_filepath = format!("{}.fasta", output_filepath.rsplit_once(".").unwrap().0);

    let shannon_entroy_output_filepath = format!("{}.shannon.txt", output_filepath);

    let pb = get_spin_pb(format!("reading {}", filepath), DEFAULT_INTERVAL);
    let all_records = reader
        .records()
        .into_iter()
        .map(|rec| {
            pb.inc(1);
            rec
        })
        .filter(|rec| rec.is_ok())
        .map(|rec| rec.unwrap())
        .filter(|rec| !rec.is_unmapped())
        .filter(|rec| !rec.is_secondary())
        .filter(|rec| !rec.is_supplementary())
        .collect::<Vec<_>>();
    pb.finish();

    let major_pos_ins = compute_max_ins_of_each_ref_position(&all_records, Some(0), None);

    // println!("major_pos_ins:{:?}", major_pos_ins);

    let mut major_pos_ins_vec = major_pos_ins
        .iter()
        .map(|(&k, &v)| (k as usize, v as usize))
        .collect::<Vec<_>>();
    major_pos_ins_vec.sort_by_key(|v| v.0);
    // (0..(ins_size+1)).into_iter().collect::<Vec<_>>())

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
        .zip(major_start_point.into_iter())
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
    records: &Vec<Record>,
    rstart: Option<usize>,
    rend: Option<usize>,
) -> HashMap<i64, i32> {
    let mut pos2ins = HashMap::new();

    let rstart = rstart.map(|v| v as i64);
    let rend = rend.map(|v| v as i64);

    for record in records {
        // let query_locus_blacklist = get_query_locus_blacklist(record, query_locus_blacklist_gen);

        let record_ext = BamRecordExt::new(record);
        let mut start = rstart.unwrap_or(record_ext.reference_start() as i64);
        let mut end = rend.unwrap_or(record_ext.reference_end() as i64);
        start = start.max(record_ext.reference_start() as i64);
        end = end.max(record_ext.reference_end() as i64);

        let mut rpos_cursor = None;
        // let mut qpos_cursor = None;
        let mut cur_ins = 0;
        let query_end = record_ext.query_alignment_end();
        let mut aligned_pair_full = record.aligned_pairs_full().collect::<Vec<_>>();
        if aligned_pair_full.len() == 0 {
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
            // if qpos.is_some() {
            //     qpos_cursor = qpos;
            // }
            if rpos_cursor.is_none() {
                continue;
            }

            if rpos_cursor.unwrap() < start {
                continue;
            }

            // print!("{},", rpos_cursor.unwrap());
            if rpos_cursor.unwrap() >= end {
                // set the last rpos max ins and then break
                let rpos_ = rpos_cursor.unwrap();
                pos2ins.entry(rpos_ - 1).or_insert(0);
                *pos2ins.get_mut(&(rpos_ - 1)).unwrap() =
                    cmp::max(*pos2ins.get(&(rpos_ - 1)).unwrap(), cur_ins);
                break;
            }

            if let Some(qpos_) = qpos {
                if qpos_ as usize >= query_end {
                    // println!("query hit end: {}", qpos_);
                    // set the last rpos max ins and then break
                    let rpos_ = rpos_cursor.unwrap();

                    pos2ins.entry(rpos_).or_insert(0);
                    *pos2ins.get_mut(&rpos_).unwrap() =
                        cmp::max(*pos2ins.get(&rpos_).unwrap(), cur_ins);
                    break;
                }
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
    fn shanon_entropy(&self) -> f32 {
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
    records: &Vec<Record>,
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

        for [qpos, rpos] in record.aligned_pairs_full() {
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
            let tt = major_pos2major_starting_point
                .get(&(rpos_cursor.unwrap() as usize))
                .unwrap()
                + minor_cursor;
            if let Some(qpos_) = qpos {
                single_record_msa[tt] = query_seq_bytes[qpos_ as usize];
            }
        }

        single_record_msa
            .iter()
            .zip(all_locus_base_counter.iter_mut())
            .for_each(|(&base, base_counter)| {
                base_counter.counter.entry(base).and_modify(|cnt| *cnt += 1);
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

    let entropy_str = all_locus_base_counter
        .iter()
        .map(|locus| locus.shanon_entropy().to_string())
        .collect::<Vec<_>>()
        .join(",");

    let mut ofile = fs::File::create(shann_entropy_ofilename).unwrap();
    ofile.write_all(entropy_str.as_bytes()).unwrap();
}
