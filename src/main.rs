use bevy::input::keyboard::KeyCode;
use bevy::input::keyboard::KeyboardInput;
use bevy::input::ButtonState;
use bevy::prelude::*;
use bevy::render::mesh::PrimitiveTopology;
use bevy_inspector_egui::WorldInspectorPlugin;
use smooth_bevy_cameras::controllers::orbit::{
  OrbitCameraBundle, OrbitCameraController, OrbitCameraPlugin,
};
use smooth_bevy_cameras::LookTransformPlugin;

use mesh_lib::vec_math::{
  cross, dot, normalize, points_cross, scale, sub, vec_length, vec_length_square,
};
use mesh_lib::{Crease, Fold};
use std::fs;

use bevy::render::mesh::{Indices, Mesh};

struct Record {
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
struct FpsText;

fn main() {
  let axial_stiffness = 20.0;

  let data = fs::read_to_string("./mesh-lib/src/bird2.fold").unwrap();
  let mut fold: Fold = serde_json::from_str(&data).unwrap();
  let creases = fold.get_creases();
  let velocity = vec![[0.0f32, 0.0f32, 0.0f32]; fold.vertices_coords.len()];

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
    .add_system(print_keyboard_event_system)
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
    material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
    ..default()
  });

  //ith multiple sections
  commands
    .spawn_bundle(
      // Create a TextBundle that has a Text with a single section.
      TextBundle::from_section(
        // Accepts a `String` or any type that converts into a `String`, such as `&str`
        "hello\nbevy!",
        TextStyle {
          font: asset_server.load("fonts/FiraSans-Bold.ttf"),
          font_size: 100.0,
          color: Color::WHITE,
        },
      ) // Set the alignment of the Text
      .with_text_alignment(TextAlignment::TOP_CENTER)
      // Set the style of the TextBundle itself.
      .with_style(Style {
        align_self: AlignSelf::FlexEnd,
        position_type: PositionType::Absolute,
        position: UiRect {
          bottom: Val::Px(5.0),
          right: Val::Px(15.0),
          ..default()
        },
        ..default()
      }),
    )
    .insert(FpsText);

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

/// This system prints out all keyboard events as they come in
fn print_keyboard_event_system(
  mut keyboard_input_events: EventReader<KeyboardInput>,
  mut record: ResMut<Record>,
) {
  for ev in keyboard_input_events.iter() {
    match ev.state {
      ButtonState::Pressed => {
        println!("Key press: {:?} ({})", ev.key_code, ev.scan_code);
        if let Some(code) = ev.key_code {
          if code == KeyCode::Space {
            if record.state == 0 {
              record.state = 1;
            } else {
              record.state = 0;
            }
          }
        }
      }
      ButtonState::Released => {}
    }
  }
}

fn joint_animation(
  mut meshes: ResMut<Assets<Mesh>>,
  mut fold_obj: ResMut<Fold>,
  mut creases: ResMut<Vec<Crease>>,
  record: Res<Record>,
  mut velocity: ResMut<Vec<[f32; 3]>>,
) {
  if record.state == 0 {
    return;
  }

  let ref_fold = &mut *fold_obj;
  let ref_creases = &mut *creases;
  let origin_face_angle = &record.face_angles;
  let edge_lengths = &record.edge_lengths;
  // calculate all normals
  let length = (ref_fold.faces_vertices).len();
  let faces_vertices = &mut ref_fold.faces_vertices;
  let positions = &mut ref_fold.vertices_coords;
  let mut f = vec![vec![0.0f32; 3]; length];

  let edges_vertices = &mut ref_fold.edges_vertices;
  for (i, idxs) in edges_vertices.iter().enumerate() {
    // edge
    let mut x01 = sub(&positions[idxs[1]], &positions[idxs[0]]);
    let new_length = vec_length(&x01);
    let k = record.axial_stiffness / edge_lengths[i];
    let force = rxn_force(k, new_length, edge_lengths[i]);
    if f32::is_nan(force) || f32::is_infinite(force) {
      panic!("err");
    }
    x01 = normalize(&x01);
    f[idxs[0]][0] -= x01[0] * force;
    f[idxs[0]][1] -= x01[1] * force;
    f[idxs[0]][2] -= x01[2] * force;

    f[idxs[1]][0] += x01[0] * force;
    f[idxs[1]][1] += x01[1] * force;
    f[idxs[1]][2] += x01[2] * force;

    let d = record.percent_damping * 2.0 * (k * 1.0).sqrt();
    let v01 = sub(&velocity[idxs[1]], &velocity[idxs[0]]);
    f[idxs[0]][0] += v01[0] * d;
    f[idxs[0]][1] += v01[1] * d;
    f[idxs[0]][2] += v01[2] * d;

    f[idxs[1]][0] -= v01[0] * d;
    f[idxs[1]][1] -= v01[1] * d;
    f[idxs[1]][2] -= v01[2] * d;

    if idxs[0] == 0 || idxs[1] == 0 {
      println!(
        "force from edge,  force={}, vec={}",
        force,
        vec_length(&[f[0][0], f[0][1], f[0][2]])
      );
    }
  }

  for (_ci, crease) in ref_creases.iter().enumerate() {
    // crease
    let vertices_idxs = crease.top_vertices_idxs;
    let theta = crease.get_theta(positions, faces_vertices);
    let mut diff = theta - record.fold_ratio * crease.target_angle;

    if vec_length(&crease.get_edge_vector(positions)) < 0.00001 {
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

    let [normal0, normal1] = crease.get_normals(positions, faces_vertices);

    let [c00, c01, h0, h1] = crease.get_0_coef(&positions);
    let [c10, c11, _h00, _h11] = crease.get_1_coef(&positions);

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
    let node0_f = scale(&normal0, rxn_force_scale / h0);
    let node1_f = scale(&normal1, rxn_force_scale / h1);

    f[vertices_idxs[0]][0] -= 1.0 * node0_f[0];
    f[vertices_idxs[0]][1] -= 1.0 * node0_f[1];
    f[vertices_idxs[0]][2] -= 1.0 * node0_f[2];

    f[vertices_idxs[1]][0] -= 1.0 * node1_f[0];
    f[vertices_idxs[1]][1] -= 1.0 * node1_f[1];
    f[vertices_idxs[1]][2] -= 1.0 * node1_f[2];

    // f[edge_vertices_idxs[0]][0] +=
    //     c10 / (c00 + c10) * node0_f[0] + c11 / (c01 + c11) * node1_f[0];
    // f[edge_vertices_idxs[0]][1] +=
    //     c10 / (c00 + c10) * node0_f[1] + c11 / (c01 + c11) * node1_f[1];
    // f[edge_vertices_idxs[0]][2] +=
    //     c10 / (c00 + c10) * node0_f[2] + c11 / (c01 + c11) * node1_f[2];

    // f[edge_vertices_idxs[1]][0] +=
    //     c00 / (c00 + c10) * node0_f[0] + c01 / (c01 + c11) * node1_f[0];
    // f[edge_vertices_idxs[1]][1] +=
    //     c00 / (c00 + c10) * node0_f[1] + c01 / (c01 + c11) * node1_f[1];
    // f[edge_vertices_idxs[1]][2] +=
    //     c00 / (c00 + c10) * node0_f[2] + c01 / (c01 + c11) * node1_f[2];

    f[edge_vertices_idxs[0]][0] += (1.0 - c00) * node0_f[0] + (1.0 - c01) * node1_f[0];
    f[edge_vertices_idxs[0]][1] += (1.0 - c00) * node0_f[1] + (1.0 - c01) * node1_f[1];
    f[edge_vertices_idxs[0]][2] += (1.0 - c00) * node0_f[2] + (1.0 - c01) * node1_f[2];

    f[edge_vertices_idxs[1]][0] += (c00) * node0_f[0] + c01 * node1_f[0];
    f[edge_vertices_idxs[1]][1] += (c00) * node0_f[1] + c01 * node1_f[1];
    f[edge_vertices_idxs[1]][2] += (c00) * node0_f[2] + c01 * node1_f[2];

    if _ci == 34 {
      println!(
        "ci={},faces={:?} ,diff={:.3}, h0={:.3},h1={:.3}",
        _ci, crease.face_idxs, diff, h0, h1
      );

      println!(
        "force from angle,  force={}",
        vec_length(&[f[0][0], f[0][1], f[0][2]])
      );
    }
    // f[edge_vertices_idxs[0]][0] += (c10) * node0_f[0] + (c11) * node1_f[0];
    // f[edge_vertices_idxs[0]][1] += (c10) * node0_f[1] + (c11) * node1_f[1];
    // f[edge_vertices_idxs[0]][2] += (c10) * node0_f[2] + (c11) * node1_f[2];

    // f[edge_vertices_idxs[1]][0] += (1.0 - c10) * node0_f[0] + (1.0 - c11) * node1_f[0];
    // f[edge_vertices_idxs[1]][1] += (1.0 - c10) * node0_f[1] + (1.0 - c11) * node1_f[1];
    // f[edge_vertices_idxs[1]][2] += (1.0 - c10) * node0_f[2] + (1.0 - c11) * node1_f[2];
  }

  for (fi, idxs) in faces_vertices.iter().enumerate() {
    // let mut ret_vec: Vec<[f64; 3]> = Vec::new();
    let a = positions[idxs[0]];
    let b = positions[idxs[1]];
    let c = positions[idxs[2]];
    let len_ab = vec_length(&sub(&a, &b));
    let len_bc = vec_length(&sub(&b, &c));
    let len_ca = vec_length(&sub(&c, &a));

    if len_ab < 0.000001 || len_bc < 0.000001 || len_ca < 0.000001 {
      continue;
    }

    let ab = normalize(&sub(&b, &a));
    let ac = normalize(&sub(&c, &a));
    let bc = normalize(&sub(&c, &b));
    let angles = [
      dot(&ab, &ac).acos(),
      (-1.0 * dot(&ab, &bc)).acos(),
      dot(&ac, &bc).acos(),
    ];

    let normal = normalize(&points_cross(&a, &b, &c));

    let diff = sub(&angles, &origin_face_angle[fi]);
    let force = scale(&diff, -1.0 * record.face_stiffness);

    let tmp_ba = scale(
      &cross(&normal, &sub(&a, &b)),
      vec_length_square(&sub(&a, &b)),
    );
    let tmp_ab = scale(&tmp_ba, -1.0);

    let tmp_bc = scale(
      &cross(&normal, &sub(&c, &b)),
      vec_length_square(&sub(&c, &b)),
    );
    let tmp_cb = scale(&tmp_bc, -1.0);

    let tmp_ca = scale(
      &cross(&normal, &sub(&a, &c)),
      vec_length_square(&sub(&a, &c)),
    );
    let tmp_ac = scale(&tmp_ca, -1.0);

    // force[1] = 0.0;
    // force[0] = 0.0;
    // force[2] = 0.0;

    f[idxs[0]][0] +=
      force[1] * tmp_ba[0] + force[0] * (tmp_ab[0] - tmp_ac[0]) - force[2] * tmp_ca[0];
    f[idxs[0]][1] +=
      force[1] * tmp_ba[1] + force[0] * (tmp_ab[1] - tmp_ac[1]) - force[2] * tmp_ca[1];
    f[idxs[0]][2] +=
      force[1] * tmp_ba[2] + force[0] * (tmp_ab[2] - tmp_ac[2]) - force[2] * tmp_ca[2];

    f[idxs[1]][0] +=
      force[1] * (tmp_bc[0] - tmp_ba[0]) + force[0] * -1.0 * tmp_ab[0] + force[2] * tmp_cb[0];
    f[idxs[1]][1] +=
      force[1] * (tmp_bc[1] - tmp_ba[1]) + force[0] * -1.0 * tmp_ab[1] + force[2] * tmp_cb[1];
    f[idxs[1]][2] +=
      force[1] * (tmp_bc[2] - tmp_ba[2]) + force[0] * -1.0 * tmp_ab[2] + force[2] * tmp_cb[2];

    f[idxs[2]][0] +=
      force[1] * -1.0 * tmp_bc[0] + force[0] * tmp_ac[0] + force[2] * (tmp_ca[0] - tmp_cb[0]);
    f[idxs[2]][1] +=
      force[1] * -1.0 * tmp_bc[1] + force[0] * tmp_ac[1] + force[2] * (tmp_ca[1] - tmp_cb[1]);
    f[idxs[2]][2] +=
      force[1] * -1.0 * tmp_bc[2] + force[0] * tmp_ac[2] + force[2] * (tmp_ca[2] - tmp_cb[2]);
  }

  //println!("{:?}", edge_lengths);

  // let edge = &fold_obj.edges_vertices;
  //let positions = &mut *fold_obj.vertices_coords;
  let delta_t = record.dt / 3.0;
  for (i, position) in &mut positions.iter_mut().enumerate() {
    //let a0 = f[i][0] / 1.0;

    if i == 0 {
      println!(
        "i={}, force={}",
        i,
        vec_length(&[f[i][0], f[i][1], f[i][2]])
      );
    }

    velocity[i][0] = velocity[i][0] + f[i][0] / 1.0 * delta_t;
    velocity[i][1] = velocity[i][1] + f[i][1] / 1.0 * delta_t;
    velocity[i][2] = velocity[i][2] + f[i][2] / 1.0 * delta_t;

    position[0] += velocity[i][0] * delta_t;
    position[1] += velocity[i][1] * delta_t;
    position[2] += velocity[i][2] * delta_t;
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
