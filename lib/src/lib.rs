
use serde::{Deserialize, Deserializer, Serialize};

#[derive(Serialize, Deserialize,Debug)]
struct Fold {
    file_spec: f32,
    file_creator: String,
    file_author: String,
    file_classes: Vec<String>,
    frame_attributes: Vec<String>,
    frame_unit: String,
    vertices_coords: Vec<[f32;3]>,
    edges_vertices:  Vec<[i32;2]>,
    edges_assignment: Vec<String>,
    faces_vertices: Vec<[i32;3]>,
}

#[cfg(test)]
mod tests {
    use std::fs;
    use super::*;

    #[test]
    fn test_parser() {
        let data = fs::read_to_string("./src/crand.fold").unwrap();
        let p: Fold = serde_json::from_str(&data).unwrap();
        println!("{:?}", p);
    }

}
