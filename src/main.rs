use bevy::app::ScheduleRunnerPlugin;
use bevy::prelude::*;
use bevy::render::mesh::PrimitiveTopology;
use bevy_inspector_egui::WorldInspectorPlugin;
use smooth_bevy_cameras::controllers::orbit::{
    OrbitCameraBundle, OrbitCameraController, OrbitCameraPlugin,
};
use smooth_bevy_cameras::LookTransformPlugin;

use mesh_lib::vec_math::{normalize, points_cross, scale, vec_length};
use mesh_lib::{Crease, Fold};
use std::fs;
use std::time::Duration;

use bevy::render::mesh::{Indices, Mesh, VertexAttributeValues};

fn main() {
    let data = fs::read_to_string("./mesh-lib/src/crand.fold").unwrap();
    let fold: Fold = serde_json::from_str(&data).unwrap();
    let creases = fold.get_creases();
    let edge_lengths = fold.get_edge_length();
    let positions = &fold.vertices_coords;
    let velocity = vec![[0.0f32, 0.0f32, 0.0f32]; positions.len()];

    App::new()
        .insert_resource(fold)
        .insert_resource(creases)
        .insert_resource(edge_lengths)
        .insert_resource(velocity)
        .add_plugins(DefaultPlugins)
        .add_plugin(LookTransformPlugin)
        .add_plugin(WorldInspectorPlugin::new())
        .add_startup_system(setup)
        // .add_stage_after(stage::UPDATE, "fixed_update", SystemStage::parallel()
        // .with_run_criteria(FixedTimestep::step(0.4))
        .add_system(joint_animation)
        .add_plugin(OrbitCameraPlugin::default())
        //.add_plugin(ScheduleRunnerPlugin(Duration::from_secs_f64(1.0 / 60.0)))
        .run();
}

/// set up a simple 3D scene
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    fold_obj: Res<Fold>,
) {
    let positions = &fold_obj.vertices_coords;
    let mut indices = Vec::with_capacity(fold_obj.faces_vertices.len() * 3);
    let normals = vec![[1.0, 1.0, 0.0]; positions.len()];
    let uvs = vec![[0.0, 0.0]; positions.len()];

    for face in &fold_obj.faces_vertices {
        for j in 0..3 {
            indices.push(face[j] as u32);
            indices.push(face[(j + 1) % 3] as u32);
        }
    }

    let mut mesh = Mesh::new(PrimitiveTopology::LineList);
    mesh.set_indices(Some(Indices::U32(indices)));
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions.clone());
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    //mesh.insert_attribute(Mesh::ATTRIBUTE_JOINT_WEIGHT, vec![1.0,1.0,1.0,1.0]);
    //mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, VertexAttributeValues::Float32x2(uvs));
    // add entities to the world

    // commands.spawn_bundle(PbrBundle {
    //     mesh: meshes.add(Mesh::from(shape::Plane { size: 5.0 })),
    //     material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
    //     ..default()
    // });

    // plane
    commands.spawn_bundle(PbrBundle {
        mesh: meshes.add(mesh),
        material: materials.add(Color::rgb(0.3, 0.5, 0.3).into()),
        ..default()
    });
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
            Vec3::new(-2.0, 5.0, 5.0),
            Vec3::new(0., 0., 0.),
        ));
}

fn joint_animation(
    mut meshes: ResMut<Assets<Mesh>>,
    mut fold_obj: ResMut<Fold>,
    mut creases: ResMut<Vec<Crease>>,
    mut edge_lengths: ResMut<Vec<f32>>,
    mut velocity: ResMut<Vec<[f32; 3]>>,
) {
    let fold_ratio = 1.0;
    let crease_stiffness = 1.9;
    let axial_stiffness = 20.0;
    let ref_fold = &mut *fold_obj;
    let ref_creases = &mut *creases;
    // calculate all normals
    let length = (ref_fold.faces_vertices).len();
    let faces_vertices = &mut ref_fold.faces_vertices;
    let positions = &mut ref_fold.vertices_coords;
    let mut normals = vec![vec![0.0f32; 3]; length];
    let mut f = vec![vec![0.0f32; 3]; length];

    // for (i, face) in faces_vertices.iter_mut().enumerate() {
    //     let a = positions.get(face[0] as usize).unwrap();
    //     let b = positions.get(face[1] as usize).unwrap();
    //     let c = positions.get(face[2] as usize).unwrap();
    //     let normal = normalize(&points_cross(&a, &b, &c));
    //     normals[i][0] = normal[0];
    //     normals[i][1] = normal[1];
    //     normals[i][2] = normal[2];
    // }

    let edges_vertices = &mut ref_fold.edges_vertices;
    for (i, idxs) in edges_vertices.iter().enumerate() {
        let v0 = positions[idxs[0]];
        let v1 = positions[idxs[1]];
        let mut v01 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];

        let new_length = vec_length(&v01);
        let diff = new_length - edge_lengths[i];

        let force = axial_stiffness * diff / edge_lengths[i];
        println!(
            "v0={:?}, v1={:?},v01={},{},{},new={}, origin={}, force={}",
            v0, v1, v01[0], v01[1], v01[2], new_length, edge_lengths[i], force
        );

        if f32::is_nan(force) || f32::is_infinite(force) {
            panic!("err");
        }
        v01 = normalize(&v01);
        f[idxs[0]][0] += v01[0] * force;
        f[idxs[0]][1] += v01[1] * force;
        f[idxs[0]][2] += v01[2] * force;

        f[idxs[1]][0] -= v01[0] * force;
        f[idxs[1]][1] -= v01[1] * force;
        f[idxs[1]][2] -= v01[2] * force;
    }

    for crease in ref_creases {
        let normals = crease.get_normals(positions, faces_vertices);
        let vertices_idxs = crease.top_vertices_idxs;
        let diff = crease.get_theta(positions, faces_vertices) - fold_ratio * crease.target_angle;
        let face_dist = crease.init_face_dist;
        let rxn_force_scale = crease.origin_lentgh * crease_stiffness * diff * 1.0;

        let node1_f = scale(
            &normals[0],
            1.0 / (vec_length(&normals[0]) / vec_length(&crease.get_edge_vector(positions)))
                * rxn_force_scale,
        );
        f[vertices_idxs[0]][0] -= node1_f[0];
        f[vertices_idxs[0]][1] -= node1_f[1];
        f[vertices_idxs[0]][2] -= node1_f[2];

        let node2_f = scale(
            &normals[1],
            1.0 / (vec_length(&normals[1]) / vec_length(&crease.get_edge_vector(positions)))
                * rxn_force_scale,
        );
        f[vertices_idxs[1]][0] -= node2_f[0];
        f[vertices_idxs[1]][1] -= node2_f[1];
        f[vertices_idxs[1]][2] -= node2_f[2];

        let edge_vertices_idxs = crease.edge_vertices_idxs;

        f[edge_vertices_idxs[0]][0] += (node1_f[0] + node2_f[0]) / 2.0;
        f[edge_vertices_idxs[0]][1] += (node1_f[1] + node2_f[1]) / 2.0;
        f[edge_vertices_idxs[0]][2] += (node1_f[2] + node2_f[2]) / 2.0;

        f[edge_vertices_idxs[1]][0] += (node1_f[0] + node2_f[0]) / 2.0;
        f[edge_vertices_idxs[1]][1] += (node1_f[1] + node2_f[1]) / 2.0;
        f[edge_vertices_idxs[1]][2] += (node1_f[2] + node2_f[2]) / 2.0;
    }

    //println!("{:?}", edge_lengths);

    // let edge = &fold_obj.edges_vertices;
    //let positions = &mut *fold_obj.vertices_coords;
    let delta_t = 1.0 / 30.0;
    let decay = 0.9;
    for (i, position) in &mut positions.iter_mut().enumerate() {
        //let a0 = f[i][0] / 1.0;

        velocity[i][0] = velocity[i][0] * decay + f[i][0] / 32.0 * delta_t;
        velocity[i][1] = velocity[i][1] * decay + f[i][1] / 32.0 * delta_t;
        velocity[i][2] = velocity[i][2] * decay + f[i][2] / 32.0 * delta_t;

        position[0] += velocity[i][0] * delta_t;
        position[1] += velocity[i][1] * delta_t;
        position[2] += velocity[i][2] * delta_t;
    }

    for (_handle_id, mesh) in meshes.iter_mut() {
        //println!("{:?}", handle_id)
        //mesh.attribute(id)
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions.to_vec());
    }
    //println!("{:?}", fold_obj)
}
