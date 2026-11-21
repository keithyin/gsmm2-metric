use std::{collections::HashMap, sync::Arc};

use crate::adapter::Mapping;
use crate::global_data::GlobalData;
use crate::kits::ReadInfo;

pub mod hp_metric_v2;
pub mod hp_tr_metric_v2;
pub mod hp_tr_metric_v3;
pub mod hp_tr_tools;
pub mod realign;

pub trait TMetric: Send + 'static + Default {
    fn csv_header() -> String;

    fn set_qname(&mut self, qname: String);
    fn get_qname(&self) -> &str;

    fn set_global_data(&mut self, global_data: Arc<GlobalData>);
    fn get_global_data(&self) -> &GlobalData;
    fn compute_metric(&mut self, read_info: &ReadInfo, reference_anchored: bool);

    fn build_metric_str(&mut self) -> String;
    fn set_target_name(&mut self, target_name: Arc<String>);
    fn get_target_name(&self) -> &str;

    fn get_mappings_mut(&mut self) -> &mut Vec<Mapping>;

    fn get_metric_str_mut(&mut self) -> &mut Option<String>;
    fn get_metric_str(&self) -> &Option<String>;

    fn set_mappings(&mut self, mappings: Vec<Mapping>) {
        *self.get_mappings_mut() = mappings;
    }

    fn add_mapping(&mut self, mapping: Mapping) {
        self.get_mappings_mut().push(mapping);
    }

    fn set_metric_str(&mut self) {
        let metric_str = self.build_metric_str();
        *self.get_metric_str_mut() = Some(metric_str);
    }
}

lazy_static::lazy_static! {
    pub static ref BASE_COMP: HashMap<u8, u8> = {
        let mut  m = HashMap::new();
        m.insert(b'A', b'T');
        m.insert(b'C', b'G');
        m.insert(b'G', b'C');
        m.insert(b'T', b'A');
        m
    };
}
