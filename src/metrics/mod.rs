use std::sync::Arc;

use crate::global_data::GlobalData;

pub mod hp_tr_metric;

pub trait TMetric: Send + 'static + Default {
    fn csv_header() -> String;

    fn set_qname(&mut self, qname: String);
    fn get_qname(&self) -> &str;

    fn set_global_data(&mut self, global_data: Arc<GlobalData>);
    fn get_global_data(&self) -> &GlobalData;
    fn compute_metric(&mut self, read_info: &mm2::gskits::ds::ReadInfo);

    fn build_metric_str(&mut self) -> String;
    fn set_target_name(&mut self, target_name: Arc<String>);
    fn get_target_name(&self) -> &str;

    fn get_mappings_mut(&mut self) -> &mut Vec<mm2::minimap2::Mapping>;

    fn get_metric_str_mut(&mut self) -> &mut Option<String>;
    fn get_metric_str(&self) -> &Option<String>;

    fn set_mappings(&mut self, mappings: Vec<mm2::minimap2::Mapping>) {
        *self.get_mappings_mut() = mappings;
    }

    fn add_mapping(&mut self, mapping: mm2::minimap2::Mapping) {
        self.get_mappings_mut().push(mapping);
    }

    fn set_metric_str(&mut self) {
        let metric_str = self.build_metric_str();
        *self.get_metric_str_mut() = Some(metric_str);
    }
}
