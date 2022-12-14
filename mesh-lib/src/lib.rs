use bevy::prelude::Vec3;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub mod vec_math;

use crate::vec_math::points_cross_vec3;

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
  pub vertices_coords: Vec<Vec3>,
  pub edges_vertices: Vec<[usize; 2]>,
  pub edges_assignment: Vec<String>,
  pub faces_vertices: Vec<[usize; 3]>,
  dt: Option<f32>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Crease {
  // edge origin point and length
  pub edge_idx: usize,
  pub edge_vertices_idxs: [usize; 2],
  pub origin_lentgh: f32,
  pub target_angle: f32,
  pub assignment: String,
  // face index
  pub face_idxs: [usize; 2],

  // top vertice idx
  pub top_vertices_idxs: [usize; 2],
}

impl Crease {
  pub fn get_0_coef(&self, vertices_coords: &Vec<Vec3>) -> [f32; 6] {
    let p0 = vertices_coords[self.edge_vertices_idxs[0]];
    let t0 = vertices_coords[self.top_vertices_idxs[0]];
    let t1 = vertices_coords[self.top_vertices_idxs[1]];

    let mut crease = self.get_edge_vector(vertices_coords);
    let crease_length = crease.length();
    crease = crease.normalize();

    let v0 = (t0 - p0).normalize();
    let v1 = (t1 - p0).normalize();

    let cos0 = v0.dot(crease).clamp(-1.0, 1.0);
    let cos1 = v1.dot(crease).clamp(-1.0, 1.0);

    let sin0 = (1.0 - cos0 * cos0).sqrt();
    let sin1 = (1.0 - cos1 * cos1).sqrt();

    [
      cos0,
      cos1,
      (t0 - p0).length() * sin0 * 1.0,
      (t1 - p0).length() * sin1 * 1.0,
      (t0 - p0).length() * cos0 / crease_length,
      (t1 - p0).length() * cos1 / crease_length,
    ]
  }

  pub fn get_1_coef(&self, vertices_coords: &Vec<Vec3>) -> [f32; 4] {
    let p0 = vertices_coords[self.edge_vertices_idxs[1]];
    let t0 = vertices_coords[self.top_vertices_idxs[0]];
    let t1 = vertices_coords[self.top_vertices_idxs[1]];

    let crease = self.get_edge_vector(vertices_coords).normalize() * -1.0;

    let v0 = (t0 - p0).normalize();
    let v1 = (t1 - p0).normalize();

    let cos0 = v0.dot(crease).clamp(-1.0, 1.0);
    let cos1 = v1.dot(crease).clamp(-1.0, 1.0);

    let sin0 = (1.0 - cos0 * cos0).sqrt();
    let sin1 = (1.0 - cos1 * cos1).sqrt();

    [
      cos0,
      cos1,
      (t0 - p0).length() * sin0 * 1.0,
      (t1 - p0).length() * sin1 * 1.0,
    ]
  }

  pub fn get_edge_vector(&self, vertices_coords: &Vec<Vec3>) -> Vec3 {
    // v01 = v1 - v0
    let edge_vertices_idxs = &self.edge_vertices_idxs;
    let v1 = vertices_coords[edge_vertices_idxs[1]];
    let v0 = vertices_coords[edge_vertices_idxs[0]];
    v1 - v0
  }

  pub fn get_theta(&self, normals: &Vec<Vec3>, vertices_coords: &Vec<Vec3>) -> f32 {
    let val = &self.face_idxs;
    let normal0 = normals[val[0]];
    let normal1 = normals[val[1]];

    let dot_normals = normal0.dot(normal1).clamp(-1.0, 1.0);
    let crease_vector = self.get_edge_vector(vertices_coords).normalize();
    // let res = dot(&cross(&normal0, &crease_vector), &normal1).atan2(dot_normals);
    // res
    normal0.cross(crease_vector).dot(normal1).atan2(dot_normals)
  }
}

impl Fold {
  // face1Ind, vertInd, face2Ind, ver2Ind, edgeInd, angle
  pub fn get_edge_length(&self) -> Vec<f32> {
    let edges_vertices = &self.edges_vertices;
    let vertices_coords = &self.vertices_coords;
    let mut ret_vec: Vec<f32> = Vec::new();
    for idxs in edges_vertices.iter() {
      let v0 = vertices_coords[idxs[0]];
      let v1 = vertices_coords[idxs[1]];
      //ret_vec.push(points_length(&v0, &v1));
      ret_vec.push((v1 - v0).length());
    }
    ret_vec
  }

  pub fn get_normals(&self) -> Vec<Vec3> {
    let faces_vertices = &self.faces_vertices;
    let vertices_coords = &self.vertices_coords;
    let mut ret_vec: Vec<Vec3> = Vec::new();
    for idxs in faces_vertices {
      let normal = points_cross_vec3(
        vertices_coords[idxs[0]],
        vertices_coords[idxs[1]],
        vertices_coords[idxs[2]],
      )
      .normalize();

      // let normal = normalize(&points_cross(
      //   &vertices_coords[idxs[0]],
      //   &vertices_coords[idxs[1]],
      //   &vertices_coords[idxs[2]],
      // ));
      ret_vec.push(normal);
    }
    ret_vec
  }

  pub fn get_face_angles(&self) -> Vec<[f32; 3]> {
    let faces_vertices = &self.faces_vertices;
    let positions = &self.vertices_coords;

    let mut ret_vec: Vec<[f32; 3]> = Vec::new();
    for idxs in faces_vertices.iter() {
      let a = positions[idxs[0]];
      let b = positions[idxs[1]];
      let c = positions[idxs[2]];

      let ab = (b - a).normalize();
      let ac = (c - a).normalize();
      let bc = (c - b).normalize();
      ret_vec.push([
        ab.dot(ac).acos(),
        -1.0 * ab.dot(bc).acos(),
        ac.dot(bc).acos(),
      ]);
    }
    ret_vec
  }

  pub fn get_dt(&mut self, axial_stiffness: f32) -> f32 {
    match self.dt {
      Some(x) => x,
      None => {
        let edges_vertices = &self.edges_vertices;
        let positions = &self.vertices_coords;
        let mut max_freq = 0.0;
        for idxs in edges_vertices.iter() {
          let x0 = positions[idxs[0]];
          let x1 = positions[idxs[1]];
          let x01 = x1 - x0;

          let length = x01.length();

          let k = axial_stiffness / length;
          let natural_freq = k.sqrt();

          if natural_freq > max_freq {
            max_freq = natural_freq;
          }
        }

        let dt = 1.0 / (2.0 * std::f32::consts::PI * max_freq) * 0.9;
        self.dt = Some(dt);
        dt
      }
    }
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
        "M" => -1.0 * std::f32::consts::PI,
        "V" => std::f32::consts::PI,
        "F" => 0.0f32,
        _ => 0.0f32,
      };

      if is_crease {
        // determin who is 0, 1
        let e0_idx = faces_vertices[val[0]]
          .iter()
          .position(|&x| x == idxs[0])
          .unwrap() as i32;

        let e1_idx = faces_vertices[val[0]]
          .iter()
          .position(|&x| x == idxs[1])
          .unwrap() as i32;

        let new_val = if e1_idx - e0_idx == 1 || e1_idx - e0_idx == -2 {
          [val[0], val[1]]
        } else {
          [val[1], val[0]]
        };

        let face0_top_idx = faces_vertices[new_val[0]]
          .iter()
          .position(|&x| x != idxs[0] && x != idxs[1])
          .unwrap();

        let face1_top_idx = faces_vertices[new_val[1]]
          .iter()
          .position(|&x| x != idxs[0] && x != idxs[1])
          .unwrap();

        zero_vec.push(Crease {
          edge_vertices_idxs: idxs.clone(),
          origin_lentgh: (vertices_coords[idxs[0]] - vertices_coords[idxs[1]]).length(),
          face_idxs: [new_val[0], new_val[1]],

          assignment: edges_assignment[i].clone(),
          edge_idx: i,
          top_vertices_idxs: [
            faces_vertices[new_val[0]][face0_top_idx],
            faces_vertices[new_val[1]][face1_top_idx],
          ],
          target_angle: angle,
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
    let data = fs::read_to_string("./src/bird.fold").unwrap();
    let p: Fold = serde_json::from_str(&data).unwrap();
    let vertices_coords = &p.vertices_coords;
    let faces_vertices = &p.faces_vertices;
    //println!("{:?}", p);
    let creases = p.get_creases();

    for crease in &creases {
      let theta1 = crease.get_theta(vertices_coords, faces_vertices);
      println!(
        "edge={:?}, tops={:?}",
        crease.edge_vertices_idxs, crease.top_vertices_idxs
      );
    }
    println!("len={}", creases.len())
  }
}
