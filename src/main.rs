mod system;

use bevy::prelude::*;
use bevy::render::mesh::PrimitiveTopology;
use bevy_inspector_egui::WorldInspectorPlugin;
use smooth_bevy_cameras::controllers::orbit::{
  OrbitCameraBundle, OrbitCameraController, OrbitCameraPlugin,
};
use smooth_bevy_cameras::LookTransformPlugin;

use mesh_lib::vec_math::{points_cross_vec3, scale, sub};
use mesh_lib::{Crease, Fold};
use std::fs;

use bevy::render::mesh::{Indices, Mesh};

pub struct Record {
  face_angles: Vec<[f32; 3]>,
  edge_lengths: Vec<f32>,
  dt: f32,
  fold_ratio: f32,
  crease_crease_stiffness: f32,
  flat_crease_stiffness: f32,
  face_stiffness: f32,
  axial_stiffness: f32,
  percent_damping: f32,
  state: i32,
}

#[derive(Component)]
pub struct RatioText;

#[derive(Component)]
pub struct StatusText;

fn main() {
  let axial_stiffness = 20.0;

  let data = fs::read_to_string("./mesh-lib/src/bird.fold").unwrap();
  let mut fold: Fold = serde_json::from_str(&data).unwrap();
  let creases = fold.get_creases();
  let velocity = vec![Vec3::new(0.0, 0.0, 0.0); fold.vertices_coords.len()];

  let record = Record {
    fold_ratio: 0.3,
    axial_stiffness,
    crease_crease_stiffness: 0.70,
    flat_crease_stiffness: 0.70,
    face_stiffness: 0.2,
    percent_damping: 0.48,
    dt: fold.get_dt(axial_stiffness),
    face_angles: fold.get_face_angles(),
    edge_lengths: fold.get_edge_length(),
    state: 1,
  };

  App::new()
    .insert_resource(fold)
    .insert_resource(creases)
    .insert_resource(velocity)
    .insert_resource(record)
    .add_plugins(DefaultPlugins)
    .add_plugin(LookTransformPlugin)
    .add_plugin(WorldInspectorPlugin::new())
    .add_startup_system(setup)
    .add_system(joint_animation)
    .add_system(crate::system::print_keyboard_event_system)
    .add_system(crate::system::text_update_system)
    .add_system(bevy::window::close_on_esc)
    .add_plugin(OrbitCameraPlugin::default())
    //.add_plugin(ScheduleRunnerPlugin(Duration::from_secs_f64(1.0 / 60.0)))
    .run();
}

fn rxn_force(k: f32, value: f32, target: f32) -> f32 {
  -1.0 * k * (value - target)
}
/// set up a simple 3D scene
fn setup(
  mut commands: Commands,
  mut meshes: ResMut<Assets<Mesh>>,
  mut materials: ResMut<Assets<StandardMaterial>>,
  asset_server: Res<AssetServer>,
  fold_obj: Res<Fold>,
) {
  let positions = &fold_obj.vertices_coords;
  let mut indices = Vec::with_capacity(fold_obj.faces_vertices.len() * 3);
  let normals = vec![[1.0, 1.0, 0.0]; positions.len()];
  let uvs = vec![[1.0, 1.0]; positions.len()];

  for face in &fold_obj.faces_vertices {
    for j in 0..3 {
      indices.push(face[j] as u32);
      indices.push(face[(j + 1) % 3] as u32);
    }
  }

  let mut mesh = Mesh::new(PrimitiveTopology::LineList);
  mesh.set_indices(Some(Indices::U32(indices)));
  mesh.insert_attribute(
    Mesh::ATTRIBUTE_POSITION,
    positions
      .iter()
      .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
      .collect::<Vec<[f32; 3]>>(),
  );
  mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
  mesh.insert_attribute(
    Mesh::ATTRIBUTE_COLOR,
    vec![[1.0, 0.0, 1.0, 1.0]; positions.len()],
  );
  // mesh.insert_attribute(Mesh::ATTRIBUTE_JOINT_WEIGHT, vec![1.0, 1.0, 1.0, 1.0]);
  mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
  // add entities to the world

  // plane
  commands.spawn_bundle(PbrBundle {
    mesh: meshes.add(mesh),
    material: materials.add(Color::rgb(0.1, 0.1, 0.1).into()),
    ..default()
  });

  commands
    .spawn_bundle(
      // Create a TextBundle that has a Text with a list of sections.
      TextBundle::from_sections([
        TextSection::new(
          "fold ratio: ",
          TextStyle {
            font: asset_server.load("fonts/FiraSans-Bold.ttf"),
            font_size: 30.0,
            color: Color::WHITE,
          },
        ),
        TextSection::from_style(TextStyle {
          font: asset_server.load("fonts/FiraMono-Medium.ttf"),
          font_size: 30.0,
          color: Color::GOLD,
        }),
      ])
      .with_style(Style {
        align_self: AlignSelf::FlexEnd,
        position_type: PositionType::Absolute,
        position: UiRect {
          bottom: Val::Px(5.0),
          left: Val::Px(15.0),
          ..default()
        },
        ..default()
      }),
    )
    .insert(RatioText);

  commands
    .spawn_bundle(
      // Create a TextBundle that has a Text with a list of sections.
      TextBundle::from_sections([
        TextSection::new(
          "status:",
          TextStyle {
            font: asset_server.load("fonts/FiraSans-Bold.ttf"),
            font_size: 30.0,
            color: Color::WHITE,
          },
        ),
        TextSection::from_style(TextStyle {
          font: asset_server.load("fonts/FiraMono-Medium.ttf"),
          font_size: 30.0,
          color: Color::GOLD,
        }),
      ])
      .with_style(Style {
        align_self: AlignSelf::FlexEnd,
        position_type: PositionType::Absolute,
        position: UiRect {
          bottom: Val::Px(30.0),
          left: Val::Px(15.0),
          ..default()
        },
        ..default()
      }),
    )
    .insert(StatusText);

  // light
  commands.spawn_bundle(PointLightBundle {
    point_light: PointLight {
      intensity: 1500.0,
      shadows_enabled: true,
      ..default()
    },
    transform: Transform::from_xyz(4.0, 8.0, 4.0),
    ..default()
  });

  commands
    .spawn_bundle(Camera3dBundle::default())
    .insert_bundle(OrbitCameraBundle::new(
      OrbitCameraController::default(),
      Vec3::new(0.0, 5.0, 0.0),
      Vec3::new(0., 0., 0.),
    ));
}

fn joint_animation(
  mut meshes: ResMut<Assets<Mesh>>,
  mut fold_obj: ResMut<Fold>,
  mut creases: ResMut<Vec<Crease>>,
  record: Res<Record>,
  mut velocity: ResMut<Vec<Vec3>>,
) {
  // stop or simulation
  if record.state == 0 {
    return;
  }
  // calculate all normals
  let normals = fold_obj.get_normals();

  let ref_fold = &mut *fold_obj;
  let ref_creases = &mut *creases;
  let origin_face_angle = &record.face_angles;
  let edge_lengths = &record.edge_lengths;
  let length = (ref_fold.faces_vertices).len();
  let faces_vertices = &mut ref_fold.faces_vertices;
  let positions = &mut ref_fold.vertices_coords;
  let mut f = vec![Vec3::new(0.0, 0.0, 0.0); length];

  let edges_vertices = &mut ref_fold.edges_vertices;
  for (i, idxs) in edges_vertices.iter().enumerate() {
    // edge
    let mut x01 = positions[idxs[1]] - positions[idxs[0]];
    let new_length = x01.length();
    let k = record.axial_stiffness / edge_lengths[i];
    let force = rxn_force(k, new_length, edge_lengths[i]);
    if f32::is_nan(force) || f32::is_infinite(force) {
      panic!("err");
    }
    x01 = x01.normalize();

    f[idxs[0]] -= x01 * force;
    f[idxs[1]] += x01 * force;

    let d = record.percent_damping * 2.0 * (k * 1.0).sqrt();
    let v01 = velocity[idxs[1]] - velocity[idxs[0]];

    f[idxs[0]] += v01 * d;
    f[idxs[1]] -= v01 * d;

    if idxs[0] == 0 || idxs[1] == 0 {
      println!("force from edge,  force={}, vec={}", force, f[0].length());
    }
  }

  for (_ci, crease) in ref_creases.iter().enumerate() {
    // crease
    let vertices_idxs = crease.top_vertices_idxs;
    let theta = crease.get_theta(&normals, &positions);
    let mut diff = theta - record.fold_ratio * crease.target_angle;

    if crease.get_edge_vector(positions).length() < 0.00001 {
      continue;
    }

    if diff < -5.0 {
      diff += std::f32::consts::PI * 2.0;
    } else if diff > 5.0 {
      diff -= std::f32::consts::PI * 2.0;
    }
    let crease_stiffness = if crease.assignment == "F" {
      record.flat_crease_stiffness
    } else {
      record.crease_crease_stiffness
    };

    let edge_length = edge_lengths[crease.edge_idx];
    let k = edge_length * crease_stiffness;
    let rxn_force_scale = -1.0 * k * diff;

    let normal0 = normals[crease.face_idxs[0]];
    let normal1 = normals[crease.face_idxs[1]];

    let [c00, c01, h0, h1] = crease.get_0_coef(&positions);
    let [c10, c11, _h0, _h1] = crease.get_1_coef(&positions);

    if (c00).abs() < 0.00001 || (c01).abs() < 0.00001 {
      continue;
    }

    if (c10).abs() < 0.01 || (c11).abs() < 0.00001 {
      continue;
    }

    if h0 < 0.00001 || h1 < 0.00001 {
      continue;
    }

    let edge_vertices_idxs = crease.edge_vertices_idxs;
    let node0_f = normal0 * rxn_force_scale / h0;
    let node1_f = normal1 * rxn_force_scale / h1;

    f[vertices_idxs[0]] -= node0_f;
    f[vertices_idxs[1]] -= node1_f;

    f[edge_vertices_idxs[0]] += (1.0 - c00) * node0_f + (1.0 - c01) * node1_f;
    f[edge_vertices_idxs[1]] += (c00) * node0_f + (c01) * node1_f;

    if _ci == 34 {
      println!(
        "ci={},faces={:?} ,diff={:.3}, h0={:.3},h1={:.3}",
        _ci, crease.face_idxs, diff, h0, h1
      );

      println!("force from angle,  force={}", f[0].length());
    }
  }

  for (fi, idxs) in faces_vertices.iter().enumerate() {
    let a = positions[idxs[0]];
    let b = positions[idxs[1]];
    let c = positions[idxs[2]];
    let len_ab = (a - b).length();
    let len_bc = (b - c).length();
    let len_ca = (c - a).length();

    if len_ab < 0.000001 || len_bc < 0.000001 || len_ca < 0.000001 {
      continue;
    }

    let ab = (b - a).normalize();
    let ac = (c - a).normalize();
    let bc = (c - b).normalize();
    let angles = [
      ab.dot(ac).acos(),
      -ab.dot(bc).acos(),
      ac.dot(bc).acos(), //dot(&ac, &bc).acos(),
    ];

    let normal = points_cross_vec3(a, b, c).normalize();

    let diff = sub(&angles, &origin_face_angle[fi]);
    let mut force = scale(&diff, -1.0 * record.face_stiffness);
    // force[0] = 0.0;
    // force[1] = 0.0;
    // force[2] = 0.0;

    let tmp_ba = normal.cross(a - b) * (a - b).length_squared();
    let tmp_bc = normal.cross(c - b) * (c - b).length_squared();
    let tmp_ca = normal.cross(a - c) * (a - c).length_squared();
    let tmp_ab = -1.0 * tmp_ba;
    let tmp_cb = -1.0 * tmp_bc;
    let tmp_ac = -1.0 * tmp_ca;

    f[idxs[0]] += force[1] * tmp_ba + force[0] * (tmp_ab - tmp_ac) - force[2] * tmp_ca;
    f[idxs[1]] += force[1] * (tmp_bc - tmp_ba) - force[0] * tmp_ab + force[2] * tmp_cb;
    f[idxs[2]] += -force[1] * tmp_bc + force[0] * tmp_ac + force[2] * (tmp_ca - tmp_cb);
  }

  //let positions = &mut *fold_obj.vertices_coords;
  let delta_t = record.dt;
  for (i, position) in &mut positions.iter_mut().enumerate() {
    if i == 0 {
      println!("i={}, f2={}", i, f[i].length());
    }
    velocity[i] += delta_t * f[i] / 1.0;
    *position += velocity[i] * delta_t;
  }

  for (_handle_id, mesh) in meshes.iter_mut() {
    mesh.insert_attribute(
      Mesh::ATTRIBUTE_POSITION,
      positions
        .iter()
        .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
        .collect::<Vec<[f32; 3]>>(),
    );
  }
}
