//! 原先来自 gskits 的通用工具(本项目实际用到的子集),为移除 rust-htslib
//! 依赖就地 vendor: BAM/FASTX 读取改用 noodles,其余原样移植。

pub mod bam_ext;
pub mod dna;
pub mod fastx;
pub mod pbar;
pub mod phreq;
pub mod range;
pub mod read_info;
pub mod regions;

pub use dna::reverse_complement;
pub use pbar::{DEFAULT_INTERVAL, get_bar_pb, get_spin_pb};
pub use phreq::phreq_list_2_quality;
pub use range::Range;
pub use read_info::ReadInfo;
