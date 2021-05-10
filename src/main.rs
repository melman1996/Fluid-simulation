use pixels::{Pixels, SurfaceTexture};
use winit::dpi::{LogicalSize};
use winit::event::{Event, VirtualKeyCode};
use winit::event_loop::{ControlFlow, EventLoop};
use winit_input_helper::WinitInputHelper;

mod fluid_cube;

fn main() {
    let SIZE = 100;
    let SCREEN_HEIGHT = SIZE * 4;
    let SCREEN_WIDTH = SIZE * 4;

    let event_loop = EventLoop::new();
    let mut input = WinitInputHelper::new();
    let window = create_window("Fluid simulation", &event_loop, SCREEN_WIDTH as f64, SCREEN_HEIGHT as f64);
    
    let surface_texture = SurfaceTexture::new(SCREEN_WIDTH as u32, SCREEN_HEIGHT as u32, &window);
    let mut pixels = Pixels::new(SIZE as u32, SIZE as u32, surface_texture).unwrap();

    let mut fluid_cube = fluid_cube::FluidCube::new(SIZE, 0.0001, 0.0000001, 0.2);

    for i in SIZE/2-10..SIZE/2+10 {
        for j in SIZE/2-10..SIZE/2+10 {
            fluid_cube.add_density(i, j, 1000.0);
        }
    }
    for i in 0..SIZE {
        for j in 0..SIZE {
            fluid_cube.add_velocity(i, j, (SIZE as f64 - (i as f64 - j as f64).abs()) / 10000.0, (SIZE as f64 - (i as f64 - j as f64).abs()) / 10000.0);
        }
    }

    event_loop.run(move | event, _, control_flow | {
        if let Event::RedrawRequested(_) = event {
            fluid_cube.draw(pixels.get_frame());

            if pixels.render().map_err(|e| println!("pixels.render() failed: {}", e)).is_err() {
                *control_flow = ControlFlow::Exit;
                return;
            }
        }
        if input.update(&event) {
            if input.key_pressed(VirtualKeyCode::Escape) || input.quit() {
                *control_flow = ControlFlow::Exit;
                return;
            }
            fluid_cube.step();
            window.request_redraw();
        }
    });
}

fn create_window(title: &str, event_loop: &EventLoop<()>, width: f64, height: f64) -> winit::window::Window {
    let window = winit::window::WindowBuilder::new()
        .with_visible(false)
        .with_title(title)
        .build(&event_loop)
        .unwrap();
    
    let size = LogicalSize::new(width, height);
    window.set_inner_size(size);
    window.set_visible(true);

    window
}