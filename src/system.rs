use crate::{RatioText, Record, StatusText};
use bevy::{
  input::{keyboard::KeyboardInput, ButtonState},
  prelude::{EventReader, KeyCode, Query, Res, ResMut, With, Without},
  text::Text,
};

pub fn text_update_system(
  record: Res<Record>,
  mut query: Query<&mut Text, (With<RatioText>, Without<StatusText>)>,
  mut query2: Query<&mut Text, With<StatusText>>,
) {
  for mut text in &mut query {
    text.sections[1].value = format!("{:.2}", record.fold_ratio);
  }

  for mut text in &mut query2 {
    text.sections[1].value = if record.state == 0 {
      "stop".to_string()
    } else {
      "simulation".to_string()
    }
  }
}

pub fn print_keyboard_event_system(
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

          if code == KeyCode::Q {
            record.fold_ratio -= 0.1;
          }

          if code == KeyCode::W {
            record.fold_ratio += 0.1;
          }
        }
      }
      ButtonState::Released => {}
    }
  }
}
