use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub mod vec_math;

use crate::vec_math::{cross, dot, normalize, points_cross, points_length, vec_length};

fn map_key(a: usize, b: usize) -> String {
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
    pub edges_vertices: Vec<[usize; 2]>,
    pub edges_assignment: Vec<String>,
    pub faces_vertices: Vec<[usize; 3]>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Crease {
    // edge origin point and length
    pub edge_idx: usize,
    pub edge_vertices_idxs: [usize; 2],
    pub origin_lentgh: f32,
    pub target_angle: f32,
    pub assignment: String,
    pub init_face_dist: [f32; 2],
    // face index
    pub face_idxs: [usize; 2],

    // top vertice idx
    pub top_vertices_idxs: [usize; 2],
}

impl Crease {
    pub fn get_edge_vector(&self, vertices_coords: &Vec<[f32; 3]>) -> [f32; 3] {
        // v0 -> v1
        let edge_vertices_idxs = &self.edge_vertices_idxs;
        let a = vertices_coords[edge_vertices_idxs[1]];
        let b = vertices_coords[edge_vertices_idxs[0]];
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    // var crease = creases[j];
    //  var normal1 = normals[crease.face1Index];
    //  var normal2 = normals[crease.face2Index];
    //  var dotNormals = normal1.dot(normal2);
    //  if (dotNormals < -1.0) dotNormals = -1.0;
    //  else if (dotNormals > 1.0) dotNormals = 1.0;

    //  var creaseVector = crease.getVector().normalize();
    //  //https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
    //  var theta = Math.atan2((normal1.clone().cross(creaseVector)).dot(normal2), dotNormals);
    pub fn get_normals(
        &self,
        vertices_coords: &Vec<[f32; 3]>,
        faces_vertices: &Vec<[usize; 3]>,
    ) -> [[f32; 3]; 2] {
        let val = &self.face_idxs;

        let normal1 = normalize(&points_cross(
            &vertices_coords[faces_vertices[val[0]][0]],
            &vertices_coords[faces_vertices[val[0]][1]],
            &vertices_coords[faces_vertices[val[0]][2]],
        ));

        let normal2 = normalize(&points_cross(
            &vertices_coords[faces_vertices[val[1]][0]],
            &vertices_coords[faces_vertices[val[1]][1]],
            &vertices_coords[faces_vertices[val[1]][2]],
        ));
        [normal1, normal2]
    }

    pub fn get_theta(
        &self,
        vertices_coords: &Vec<[f32; 3]>,
        faces_vertices: &Vec<[usize; 3]>,
    ) -> f32 {
        let [normal1, normal2] = self.get_normals(vertices_coords, faces_vertices);

        let dot_normals = dot(&normal1, &normal2);
        let crease_vector = self.get_edge_vector(vertices_coords);
        dot(&cross(&normal1, &crease_vector), &normal2).atan2(dot_normals)
    }
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
    pub fn get_edge_length(&self) -> Vec<f32> {
        let edges_vertices = &self.edges_vertices;
        let vertices_coords = &self.vertices_coords;
        let mut ret_vec: Vec<f32> = Vec::new();
        for idxs in edges_vertices.iter() {
            let v0 = vertices_coords[idxs[0]];
            let v1 = vertices_coords[idxs[1]];
            ret_vec.push(points_length(&v0, &v1));
        }
        ret_vec
    }

    pub fn get_creases(&self) -> Vec<Crease> {
        let mut edge_map = HashMap::new();
        let faces_vertices = &self.faces_vertices;
        let edges_vertices = &self.edges_vertices;
        let edges_assignment = &self.edges_assignment;
        let vertices_coords = &self.vertices_coords;

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
                let face1_top_idx = faces_vertices[val[0]]
                    .iter()
                    .position(|&x| x != idxs[0] && x != idxs[1])
                    .unwrap();

                let face2_top_idx = faces_vertices[val[1]]
                    .iter()
                    .position(|&x| x != idxs[0] && x != idxs[1])
                    .unwrap();

                zero_vec.push(Crease {
                    edge_vertices_idxs: idxs.clone(),
                    origin_lentgh: points_length(
                        &vertices_coords[idxs[0]],
                        &vertices_coords[idxs[1]],
                    ),
                    face_idxs: [val[0], val[1]],

                    init_face_dist: [
                        vec_length(&points_cross(
                            &vertices_coords[faces_vertices[val[0]][0]],
                            &vertices_coords[faces_vertices[val[0]][1]],
                            &vertices_coords[faces_vertices[val[0]][2]],
                        )) / 2.0,
                        vec_length(&points_cross(
                            &vertices_coords[faces_vertices[val[1]][0]],
                            &vertices_coords[faces_vertices[val[1]][1]],
                            &vertices_coords[faces_vertices[val[1]][2]],
                        )) / 2.0,
                    ],
                    assignment: edges_assignment[i].clone(),
                    edge_idx: i,
                    top_vertices_idxs: [
                        faces_vertices[val[0]][face1_top_idx],
                        faces_vertices[val[1]][face2_top_idx],
                    ],
                    target_angle: (angle * std::f64::consts::PI / 180.0) as f32,
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
        let vertices_coords = &p.vertices_coords;
        let faces_vertices = &p.faces_vertices;
        //println!("{:?}", p);
        let creases = p.get_creases();

        for crease in creases {
            let theta1 = crease.get_theta(vertices_coords, faces_vertices);
            println!("{:?}", theta1);
        }
    }
}
