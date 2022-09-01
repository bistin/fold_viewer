use bevy::prelude::*;
use bevy::render::mesh::PrimitiveTopology;
use bevy_inspector_egui::WorldInspectorPlugin;
use smooth_bevy_cameras::controllers::orbit::{
    OrbitCameraBundle, OrbitCameraController, OrbitCameraPlugin,
};
use smooth_bevy_cameras::LookTransformPlugin;

use mesh_lib::vec_math::{normalize, points_cross};
use mesh_lib::{Crease, Fold};
use std::fs;

use bevy::render::mesh::{Indices, Mesh, VertexAttributeValues};

fn main() {
    let data = fs::read_to_string("./mesh-lib/src/crand.fold").unwrap();
    let fold: Fold = serde_json::from_str(&data).unwrap();
    let creases = fold.get_creases();
    let edge_lengths = fold.get_edge_length();

    App::new()
        .insert_resource(fold)
        .insert_resource(creases)
        .insert_resource(edge_lengths)
        .add_plugins(DefaultPlugins)
        .add_plugin(LookTransformPlugin)
        .add_plugin(WorldInspectorPlugin::new())
        .add_startup_system(setup)
        .add_system(joint_animation)
        .add_plugin(OrbitCameraPlugin::default())
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
) {
    let ref_fold = &mut *fold_obj;
    let ref_creases = &mut *creases;
    // calculate all normals
    let length = (ref_fold.faces_vertices).len();
    let faces_vertices = &mut ref_fold.faces_vertices;
    let positions = &mut ref_fold.vertices_coords;
    let mut normals = vec![vec![0.0f32; 3]; length];

    for (i, face) in faces_vertices.iter_mut().enumerate() {
        let a = positions.get(face[0] as usize).unwrap();
        let b = positions.get(face[1] as usize).unwrap();
        let c = positions.get(face[2] as usize).unwrap();
        let normal = normalize(&points_cross(&a, &b, &c));
        normals[i][0] = normal[0];
        normals[i][1] = normal[1];
        normals[i][2] = normal[2];
    }

    println!("{:?}", edge_lengths);

    let edge = &fold_obj.edges_vertices;
    let positions = &mut *fold_obj.vertices_coords;
    for i in &mut positions.iter_mut() {
        i[1] = i[1] + 0.005
    }

    for (_handle_id, mesh) in meshes.iter_mut() {
        //println!("{:?}", handle_id)
        //mesh.attribute(id)
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions.to_vec());
    }
    //println!("{:?}", fold_obj)
}
