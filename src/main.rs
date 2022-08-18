use bevy::prelude::*;
use bevy::render::mesh::{self, PrimitiveTopology};
use bevy_inspector_egui::WorldInspectorPlugin;
use mesh_lib::Fold;
use std::fs;

use bevy::{
    asset::{AddAsset, AssetLoader, LoadContext, LoadedAsset},
    prelude::*,
    render::{
        mesh::{Indices, Mesh, VertexAttributeValues},
    },
    utils::BoxedFuture,
};


fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(WorldInspectorPlugin::new())
        .add_startup_system(setup)
        .run();
}


// fn stl_to_wireframe_mesh(stl: &stl_io::IndexedMesh) -> Mesh {
//     let mut mesh = Mesh::new(PrimitiveTopology::LineList);

//     let positions = stl.vertices.iter().map(|v| [v[0], v[1], v[2]]).collect();
//     let mut indices = Vec::with_capacity(stl.faces.len() * 3);
//     let normals = vec![[1.0, 0.0, 0.0]; stl.vertices.len()];
//     let uvs = vec![[0.0, 0.0]; stl.vertices.len()];

//     for face in &stl.faces {
//         for j in 0..3 {
//             indices.push(face.vertices[j] as u32);
//             indices.push(face.vertices[(j + 1) % 3] as u32);
//         }
//     }
//     print!("{:?}", indices);
//     mesh.insert_attribute(
//         Mesh::ATTRIBUTE_POSITION,
//         VertexAttributeValues::Float32x3(positions),
//     );
//     mesh.insert_attribute(
//         Mesh::ATTRIBUTE_NORMAL,
//         VertexAttributeValues::Float32x3(normals),
//     );
//     mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, VertexAttributeValues::Float32x2(uvs));
//     mesh.set_indices(Some(Indices::U32(indices)));

//     mesh
// }

/// set up a simple 3D scene
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    ) {
   

    let data = fs::read_to_string("./mesh-lib/src/crand.fold").unwrap();
    let fold: Fold = serde_json::from_str(&data).unwrap();

    println!("{:?}", fold);

    let positions = fold.vertices_coords;
    let mut indices = Vec::with_capacity(fold.faces_vertices.len() * 3);
    let normals = vec![[0.0, 1.0, 0.0]; positions.len()];
    let uvs = vec![[0.0, 0.0]; positions.len()];

    for face in &fold.faces_vertices {
        for j in 0..3 {
            indices.push(face[j] as u32);
            indices.push(face[(j + 1) % 3] as u32);
        }
    }
    print!("{:?}", indices);




    let mut mesh = Mesh::new(PrimitiveTopology::LineList);
    mesh.set_indices(Some(Indices::U32(indices)));
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
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
    // camera
    // commands.spawn_bundle(Camera3dBundle {
    //     transform: Transform::from_xyz(0.0, 0.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    //     ..default()
    // });
    commands.spawn_bundle(Camera3dBundle {
        transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
    // commands.spawn(Camera3dBundle {
    //     transform: Transform::from_translation(Vec3::new(-2.0, 2.5, 5.0))
    //         .looking_at(Vec3::default(), Vec3::unit_y()),
    //     ..Default::default()
    // });
    // commands.spawn_bundle(PerspectiveCameraBundle {
    //     transform: Transform::from_xyz(-2.0, 2.5, 5.0)
    //         .looking_at(Vec3::ZERO, Vec3::Y),
    //         ..default()
    // });
}