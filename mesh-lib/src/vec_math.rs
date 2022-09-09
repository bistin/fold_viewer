pub fn points_cross(a: &[f32; 3], b: &[f32; 3], c: &[f32; 3]) -> [f32; 3] {
    let ab = sub(&b, &a);
    let ac = sub(&c, &a);

    cross(&ac, &ab)
}

pub fn cross(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn dot(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn sub(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn add(a: &[f32; 3], b: &[f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn normalize(vector: &[f32; 3]) -> [f32; 3] {
    let length_sqr = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    let length = length_sqr.sqrt();

    [vector[0] / length, vector[1] / length, vector[2] / length]
}

pub fn vec_length(vector: &[f32; 3]) -> f32 {
    let length_sqr = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    length_sqr.sqrt()
}

pub fn points_length(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    let tmp_vec = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    vec_length(&tmp_vec)
}

pub fn length(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    let ab = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

    (ab[0].powi(2) + ab[1].powi(2) + ab[2].powi(2)).sqrt()
}

pub fn scale(a: &[f32; 3], scalar: f32) -> [f32; 3] {
    [a[0] * scalar, a[1] * scalar, a[2] * scalar]
}
