use bevy::math::Vec3;

pub fn points_cross_vec3(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
  let ab = b - a;
  let ac = c - a;
  ac.cross(ab)
}

pub fn sub(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] {
  [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale(a: &[f32; 3], scalar: f32) -> [f32; 3] {
  [a[0] * scalar, a[1] * scalar, a[2] * scalar]
}
