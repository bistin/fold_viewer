use crate::{RatioText, Record, StatusText};
use bevy::prelude::*;

pub fn text_update_system(
  record: Res<Record>,
  mut query: Query<&Children, (With<RatioText>, Without<StatusText>)>,
  mut query2: Query<&Children, With<StatusText>>,
  mut text_query: Query<&mut TextSpan>,
) {
  for children in &mut query {
    if let Some(child) = children.first() {
      if let Ok(mut span) = text_query.get_mut(*child) {
        **span = format!("{:.2}", record.fold_ratio);
      }
    }
  }

  for children in &mut query2 {
    if let Some(child) = children.first() {
      if let Ok(mut span) = text_query.get_mut(*child) {
        **span = if record.state == 0 {
          "stop".to_string()
        } else {
          "simulation".to_string()
        };
      }
    }
  }
}

pub fn print_keyboard_event_system(
  keyboard_input: Res<ButtonInput<KeyCode>>,
  mut record: ResMut<Record>,
) {
  if keyboard_input.just_pressed(KeyCode::Space) {
    if record.state == 0 {
      record.state = 1;
    } else {
      record.state = 0;
    }
  }

  if keyboard_input.just_pressed(KeyCode::KeyQ) {
    record.fold_ratio -= 0.1;
  }

  if keyboard_input.just_pressed(KeyCode::KeyW) {
    record.fold_ratio += 0.1;
  }
}
