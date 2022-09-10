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
    face_angles: Vec<[f64; 3]>,
    dt: f64,
}

fn main() {
    let axial_stiffness = 20.0;
    let data = fs::read_to_string("./mesh-lib/src/bird.fold").unwrap();
    let mut fold: Fold = serde_json::from_str(&data).unwrap();
    let creases = fold.get_creases();
    let edge_lengths = fold.get_edge_length();
    // init face_angles
    let face_angles = fold.get_face_angles();
    // let positions = &fold.vertices_coords;
    let velocity = vec![[0.0f64, 0.0f64, 0.0f64]; fold.vertices_coords.len()];
    let dt = fold.get_dt(axial_stiffness);

    App::new()
        .insert_resource(fold)
        .insert_resource(creases)
        .insert_resource(edge_lengths)
        .insert_resource(velocity)
        .insert_resource(Record { face_angles, dt })
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

fn rxn_force(k: f64, value: f64, target: f64) -> f64 {
    -1.0 * k * (value - target)
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
    mesh.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        positions
            .iter()
            .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
            .collect::<Vec<[f32; 3]>>(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    // mesh.insert_attribute(Mesh::ATTRIBUTE_JOINT_WEIGHT, vec![1.0, 1.0, 1.0, 1.0]);
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, uvs);
    // add entities to the world

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
    record: Res<Record>,
    edge_lengths: Res<Vec<f64>>,
    mut velocity: ResMut<Vec<[f64; 3]>>,
) {
    let fold_ratio = 1.0;
    let crease_crease_stiffness = 0.7;
    let flat_crease_stiffness = 0.9;
    let face_stiffness = 0.2;
    let axial_stiffness = 20.0;
    let percent_damping = 0.45;
    let ref_fold = &mut *fold_obj;
    let ref_creases = &mut *creases;
    let origin_face_angle = &record.face_angles;
    // calculate all normals
    let length = (ref_fold.faces_vertices).len();
    let faces_vertices = &mut ref_fold.faces_vertices;
    let positions = &mut ref_fold.vertices_coords;
    let mut f = vec![vec![0.0f64; 3]; length];

    let edges_vertices = &mut ref_fold.edges_vertices;
    for (i, idxs) in edges_vertices.iter().enumerate() {
        // edge
        let mut x01 = sub(&positions[idxs[1]], &positions[idxs[0]]);
        let new_length = vec_length(&x01);
        let k = axial_stiffness / edge_lengths[i];
        let force = rxn_force(k, new_length, edge_lengths[i]);
        if f64::is_nan(force) || f64::is_infinite(force) {
            panic!("err");
        }
        x01 = normalize(&x01);
        f[idxs[0]][0] -= x01[0] * force;
        f[idxs[0]][1] -= x01[1] * force;
        f[idxs[0]][2] -= x01[2] * force;

        f[idxs[1]][0] += x01[0] * force;
        f[idxs[1]][1] += x01[1] * force;
        f[idxs[1]][2] += x01[2] * force;

        let d = percent_damping * 2.0 * (k * 1.0).sqrt();
        let v01 = sub(&velocity[idxs[1]], &velocity[idxs[0]]);
        f[idxs[0]][0] += v01[0] * d;
        f[idxs[0]][1] += v01[1] * d;
        f[idxs[0]][2] += v01[2] * d;

        f[idxs[1]][0] -= v01[0] * d;
        f[idxs[1]][1] -= v01[1] * d;
        f[idxs[1]][2] -= v01[2] * d;
    }

    for (_ci, crease) in ref_creases.iter().enumerate() {
        // crease
        let vertices_idxs = crease.top_vertices_idxs;
        let theta = crease.get_theta(positions, faces_vertices);
        let mut diff = theta - fold_ratio * crease.target_angle;

        if diff.abs() < 0.0001 {
            diff = 0.0;
        }

        if diff < -5.0 {
            diff += std::f64::consts::PI * 2.0;
        } else if diff > 5.0 {
            diff -= std::f64::consts::PI * 2.0;
        }
        // theta = diff + fold_ratio * crease.target_angle;
        let crease_stiffness = if crease.assignment == "F" {
            flat_crease_stiffness
        } else {
            crease_crease_stiffness
        };

        let edge_length = edge_lengths[crease.edge_idx];
        let k = edge_length * crease_stiffness;
        let rxn_force_scale = -1.0 * k * diff;

        let [normal0, normal1] = crease.get_normals(positions, faces_vertices);

        let [c00, c01, h0, h1] = crease.get_0_coef(&positions);
        let [c10, c11, _h00, _h11] = crease.get_1_coef(&positions);

        if _ci > 0 {
            println!(
                "theta={}, target = {}, diff={}, ci={}, force={}, h1={}, h2={}",
                theta,
                fold_ratio * crease.target_angle,
                diff,
                _ci,
                rxn_force_scale,
                h0,
                h1
            );
        }

        // if h1 < 0.00001 || h0 < 0.00001 {
        //     continue;
        // }

        let edge_vertices_idxs = crease.edge_vertices_idxs;
        let node0_f = scale(&normal0, 1.0 / (h0) * rxn_force_scale);
        let node1_f = scale(&normal1, 1.0 / (h1) * rxn_force_scale);

        f[vertices_idxs[0]][0] -= node0_f[0];
        f[vertices_idxs[0]][1] -= node0_f[1];
        f[vertices_idxs[0]][2] -= node0_f[2];

        f[vertices_idxs[1]][0] -= node1_f[0];
        f[vertices_idxs[1]][1] -= node1_f[1];
        f[vertices_idxs[1]][2] -= node1_f[2];

        f[edge_vertices_idxs[0]][0] +=
            c10 / (c00 + c10) * node0_f[0] + c11 / (c01 + c11) * node1_f[0];
        f[edge_vertices_idxs[0]][1] +=
            c10 / (c00 + c10) * node0_f[1] + c11 / (c01 + c11) * node1_f[1];
        f[edge_vertices_idxs[0]][2] +=
            c10 / (c00 + c10) * node0_f[2] + c11 / (c01 + c11) * node1_f[2];

        f[edge_vertices_idxs[1]][0] +=
            c00 / (c00 + c10) * node0_f[0] + c01 / (c01 + c11) * node1_f[0];
        f[edge_vertices_idxs[1]][1] +=
            c00 / (c00 + c10) * node0_f[1] + c01 / (c01 + c11) * node1_f[1];
        f[edge_vertices_idxs[1]][2] +=
            c00 / (c00 + c10) * node0_f[2] + c01 / (c01 + c11) * node1_f[2];
    }

    for (fi, idxs) in faces_vertices.iter().enumerate() {
        // let mut ret_vec: Vec<[f64; 3]> = Vec::new();
        let a = positions[idxs[0]];
        let b = positions[idxs[1]];
        let c = positions[idxs[2]];
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
        let mut force = scale(&diff, -1.0 * face_stiffness);

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

        //force[1] = 0.0;
        //force[0] = 0.0;
        //force[2] = 0.0;

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

        // nominalTriangles[4 * i] = Math.acos(ab.dot(ac));
        // nominalTriangles[4 * i + 1] = Math.acos(-1 * ab.dot(bc));
        // nominalTriangles[4 * i + 2] = Math.acos(ac.dot(bc));
    }

    //println!("{:?}", edge_lengths);

    // let edge = &fold_obj.edges_vertices;
    //let positions = &mut *fold_obj.vertices_coords;
    let delta_t = record.dt;
    let decay = 1.0;
    for (i, position) in &mut positions.iter_mut().enumerate() {
        //let a0 = f[i][0] / 1.0;

        velocity[i][0] = velocity[i][0] * decay + f[i][0] / 1.0 * delta_t;
        velocity[i][1] = velocity[i][1] * decay + f[i][1] / 1.0 * delta_t;
        velocity[i][2] = velocity[i][2] * decay + f[i][2] / 1.0 * delta_t;

        position[0] += velocity[i][0] * delta_t;
        position[1] += velocity[i][1] * delta_t;
        position[2] += velocity[i][2] * delta_t;
    }

    for (_handle_id, mesh) in meshes.iter_mut() {
        //println!("{:?}", handle_id)
        //mesh.attribute(id)
        mesh.insert_attribute(
            Mesh::ATTRIBUTE_POSITION,
            positions
                .iter()
                .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
                .collect::<Vec<[f32; 3]>>(),
        );
    }
    //println!("{:?}", fold_obj)
}
