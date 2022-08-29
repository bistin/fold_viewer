use serde::{Deserialize, Serialize};
use std::collections::HashMap;

fn map_key(a: i32, b: i32) -> String {
    if a < b {
        format!("{}_{}", a, b)
    } else {
        format!("{}_{}", b, a)
    }
}

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

#[derive(Serialize, Deserialize, Debug)]
pub struct Crease {
    pub vertices_idxs: [i32; 2],
    pub top_idxs: [i32; 2],
    pub assignment: String,
    pub edge_idx: i32,
    pub target_angle: f64,
}

// for (var i=0;i<fold.edges_assignment.length;i++){
//      var assignment = fold.edges_assignment[i];
//      if (assignment == "M") foldAngles.push(-180);
//      else if (assignment == "V") foldAngles.push(180);
//      else if (assignment == "F") foldAngles.push(0);
//      else foldAngles.push(null);
//  }
impl Fold {
    // face1Ind, vertInd, face2Ind, ver2Ind, edgeInd, angle
    pub fn get_creases(&self) -> Vec<Crease> {
        let mut edge_map = HashMap::new();
        let faces_vertices = &self.faces_vertices;
        let edges_vertices = &self.edges_vertices;
        let edges_assignment = &self.edges_assignment;

        let mut zero_vec: Vec<Crease> = Vec::new();
        for (i, idxs) in faces_vertices.iter().enumerate() {
            edge_map
                .entry(map_key(idxs[0], idxs[1]))
                .or_insert(Vec::new())
                .push(i);
            edge_map
                .entry(map_key(idxs[1], idxs[2]))
                .or_insert(Vec::new())
                .push(i);

            edge_map
                .entry(map_key(idxs[0], idxs[2]))
                .or_insert(Vec::new())
                .push(i);
        }

        for (i, idxs) in edges_vertices.iter().enumerate() {
            let key = map_key(idxs[0], idxs[1]);
            let val = edge_map.get(&key).unwrap();
            let assignment = &edges_assignment[i];
            let is_crease = assignment == "M" || assignment == "V" || assignment == "F";
            let angle = match assignment.as_str() {
                "M" => -180.0f64,
                "V" => 180.0f64,
                "F" => 0.0f64,
                _ => 0.0f64,
            };

            if is_crease {
                let face1_idx = faces_vertices[val[0]]
                    .iter()
                    .position(|&x| x != idxs[0] && x != idxs[1])
                    .unwrap();

                let face2_idx = faces_vertices[val[1]]
                    .iter()
                    .position(|&x| x != idxs[0] && x != idxs[1])
                    .unwrap();

                zero_vec.push(Crease {
                    vertices_idxs: [val[0] as i32, val[1] as i32],
                    assignment: edges_assignment[i].clone(),
                    edge_idx: i as i32,
                    top_idxs: [
                        faces_vertices[val[0]][face1_idx],
                        faces_vertices[val[1]][face2_idx],
                    ],
                    target_angle: angle * std::f64::consts::PI / 180.0,
                });
            }
        }
        return zero_vec;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_parser() {
        let data = fs::read_to_string("./src/crand.fold").unwrap();
        let p: Fold = serde_json::from_str(&data).unwrap();
        //println!("{:?}", p);
        p.get_creases();
    }
}
