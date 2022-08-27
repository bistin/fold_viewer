use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Fold {
    pub file_spec: f32,
    pub file_creator: String,
    pub file_author: String,
    pub file_classes: Vec<String>,
    pub frame_attributes: Vec<String>,
    pub frame_unit: String,
    pub vertices_coords: Vec<[f32; 3]>,
    pub edges_vertices: Vec<[i32; 2]>,
    pub edges_assignment: Vec<String>,
    pub faces_vertices: Vec<[i32; 3]>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_parser() {
        let data = fs::read_to_string("./src/crand.fold").unwrap();
        let p: Fold = serde_json::from_str(&data).unwrap();
        println!("{:?}", p);
    }
}
