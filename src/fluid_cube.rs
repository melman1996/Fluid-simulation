pub struct FluidCube {
    size: usize,
    dt: f64,
    diff: f64,
    visc: f64,

    s: Vec<f64>,
    density: Vec<f64>,

    vx: Vec<f64>,
    vy: Vec<f64>,

    vx0: Vec<f64>,
    vy0: Vec<f64>,
}

impl FluidCube {
    pub fn new(size: usize, diffusion: f64, viscosity: f64, dt: f64) -> Self {
        FluidCube {
            size,
            dt,
            diff: diffusion,
            visc: viscosity,
            s: vec![0.0; size * size],
            density: vec![0.0; size * size],
            vx: vec![0.0; size * size],
            vy: vec![0.0; size * size],
            vx0: vec![0.0; size * size],
            vy0: vec![0.0; size * size],
        }
    }

    pub fn add_density(&mut self, x: usize, y: usize, amount: f64) {
        let index = IX(x, y, self.size);
        self.density[index] += amount;
    }

    pub fn add_velocity(&mut self, x: usize, y: usize, amount_x: f64, amount_y: f64) {
        let index = IX(x, y, self.size);
        self.vx[index] += amount_x;
        self.vy[index] += amount_y;
    }

    pub fn draw(&self, screen: &mut [u8]) {
        for (c, pix) in self.density.iter().zip(screen.chunks_exact_mut(4)) {
            let color = [*c as u8, *c as u8, *c as u8, *c as u8];
            pix.copy_from_slice(&color);
        }
    }

    pub fn step(&mut self){
        Self::diffuse(1, &mut self.vx0, &mut self.vx, self.visc, self.dt, 16, self.size);
        Self::diffuse(1, &mut self.vy0, &mut self.vy, self.visc, self.dt, 16, self.size);

        Self::project(&mut self.vx0, &mut self.vy0, &mut self.vx, &mut self.vy, self.size);

        Self::advect_x(1, &mut self.vx, &mut self.vx0, &mut self.vy0, self.dt, self.size);
        Self::advect_y(1, &mut self.vy, &mut self.vx0, &mut self.vy0, self.dt, self.size);

        Self::project(&mut self.vx, &mut self.vy, &mut self.vx0, &mut self.vy0, self.size);
        Self::diffuse(0, &mut self.s, &mut self.density, self.diff, self.dt, 16, self.size);
        Self::advect(0, &mut self.density, &mut self.s, &mut self.vx, &mut self.vy, self.dt, self.size);
    }

    fn diffuse(b: i32, x: &mut Vec<f64>, x0: &Vec<f64>, diff: f64, dt: f64, iter: i32, size: usize) {
        let a = dt * diff * ((size - 2) * (size - 2)) as f64;
        Self::lin_solve(b, x, x0, a, 1.0 + 4.0 * a, iter, size);
    }

    fn lin_solve(b: i32, x: &mut Vec<f64>, x0: &Vec<f64>, a: f64, c: f64, iter: i32, size: usize) {
        let cRecip = 1.0 / c;
        for t in 0..iter {
            for j in 1..size-1 {
                for i in 1..size-1 {
                    x[IX(i, j, size)] = (x0[IX(i, j, size)] + a * (x[IX(i + 1, j, size)] + x[IX(i - 1, j, size)] + x[IX(i, j + 1, size)] + x[IX(i, j - 1, size)])) * cRecip;
                }
            }
            Self::set_bnd(b, x, size);
        }
    }

    fn set_bnd(b: i32, x: &mut Vec<f64>, size: usize) {
        for i in 1..size-1 {
            x[IX(i, 0, size)] = if b == 2 {
                -x[IX(i, 1, size)]
            } else {
                x[IX(i, 1, size)]
            };
            x[IX(i, size - 1, size)] = if b == 2 {
                -x[IX(i, size - 2, size)]
            } else {
                x[IX(i, size - 2, size)]
            };
        }
        for j in 1..size-1 {
            x[IX(0, j, size)] = if b == 1 {
                -x[IX(1, j, size)]
            } else {
                x[IX(1, j, size)]
            };
            x[IX(size - 1, j, size)] = if b == 1 {
                -x[IX(size - 2, j, size)]
            } else {
                x[IX(size - 2, j, size)]
            };
        }
        x[IX(0, 0, size)] = 0.5 * (x[IX(1, 0, size)] + x[IX(0, 1, size)]);
        x[IX(0, size - 1, size)] = 0.5 * (x[IX(1, size - 1, size)] + x[IX(0, size - 2, size)]);
        x[IX(size - 1, 0, size)] = 0.5 * (x[IX(size - 2, 0, size)] + x[IX(size - 1, 1, size)]);
        x[IX(size - 1, size - 1, size)] = 0.5 * (x[IX(size - 2, size - 1, size)] + x[IX(size - 1, size - 2, size)]);
    }

    fn project(veloc_x: &mut Vec<f64>, veloc_y: &mut Vec<f64>, p: &mut Vec<f64>, div: &mut Vec<f64>, size: usize) {
        for j in 1..size-1 {
            for i in 1..size-1 {
                div[IX(i, j, size)] = (-0.5 * (veloc_x[IX(i + 1, j, size)] - veloc_x[IX(i - 1, j, size)] + veloc_y[IX(i, j + 1, size)] - veloc_y[IX(i, j - 1, size)])) / size as f64;
                p[IX(i, j, size)] = 0.0;
            }
        }
        Self::set_bnd(0, div, size);
        Self::set_bnd(0, p, size);
        Self::lin_solve(0, p, div, 1.0, 4.0, 4, size);

        for j in 1..size-1 {
            for i in 1..size-1 {
                veloc_x[IX(i, j, size)] -= 0.5 * (p[IX(i + 1, j, size)] - p[IX(i - 1, j, size)]) * size as f64;
                veloc_y[IX(i, j, size)] -= 0.5 * (p[IX(i, j + 1, size)] - p[IX(i, j - 1, size)]) * size as f64;
            }
        }

        Self::set_bnd(1, veloc_x, size);
        Self::set_bnd(2, veloc_y, size);
    }

    fn advect(
        b: i32,
        d: &mut Vec<f64>,
        d0: &mut Vec<f64>,
        veloc_x: &mut Vec<f64>,
        veloc_y: &mut Vec<f64>,
        dt: f64,
        size: usize,
    ) {
        let dt0 = dt * size as f64;
        for i in 1..size-1 {
            for j in 1..size-1 {
                let x = i as f64 - dt0 * veloc_x[IX(i, j, size)];
                let x = if x < 0.5 {
                    0.5
                } else if x > size as f64 + 0.5 {
                    size as f64 + 0.5
                } 
                else {
                    x
                };

                let y = j as f64 - dt0 * veloc_y[IX(i, j, size)];
                let y = if y < 0.5 {
                    0.5
                } else if y > size as f64 + 0.5 {
                    size as f64 + 0.5
                } 
                else {
                    y
                };

                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f64;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f64;
                let t0 = 1.0 - t1;

                d[IX(i, j, size)] = s0 * (t0 * d0[IX(i0, j0, size)] + t1 * d0[IX(i0, j1, size)])
                                      + s1 * (t0 * d0[IX(i1, j0, size)] + t1 * d0[IX(i1, j1, size)]);
            }
        }
        Self::set_bnd(b, d, size);
    }

    fn advect_x(
        b: i32,
        d: &mut Vec<f64>,
        veloc_x: &mut Vec<f64>,
        veloc_y: &mut Vec<f64>,
        dt: f64,
        size: usize,
    ) {
        let dt0 = dt * size as f64;
        for i in 1..size-1 {
            for j in 1..size-1 {
                let x = i as f64 - dt0 * veloc_x[IX(i, j, size)];
                let x = if x < 0.5 {
                    0.5
                } else if x > size as f64 + 0.5 {
                    size as f64 + 0.5
                } 
                else {
                    x
                };

                let y = j as f64 - dt0 * veloc_y[IX(i, j, size)];
                let y = if y < 0.5 {
                    0.5
                } else if y > size as f64 + 0.5 {
                    size as f64 + 0.5
                } 
                else {
                    y
                };

                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f64;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f64;
                let t0 = 1.0 - t1;

                d[IX(i, j, size)] = s0 * (t0 * veloc_x[IX(i0, j0, size)] + t1 * veloc_x[IX(i0, j1, size)])
                                      + s1 * (t0 * veloc_x[IX(i1, j0, size)] + t1 * veloc_x[IX(i1, j1, size)]);
            }
        }
        Self::set_bnd(b, d, size);
    }

    fn advect_y(
        b: i32,
        d: &mut Vec<f64>,
        veloc_x: &mut Vec<f64>,
        veloc_y: &mut Vec<f64>,
        dt: f64,
        size: usize,
    ) {
        let dt0 = dt * size as f64;
        for i in 1..size-1 {
            for j in 1..size-1 {
                let x = i as f64 - dt0 * veloc_x[IX(i, j, size)];
                let x = if x < 0.5 {
                    0.5
                } else if x > size as f64 + 0.5 {
                    size as f64 + 0.5
                } 
                else {
                    x
                };

                let y = j as f64 - dt0 * veloc_y[IX(i, j, size)];
                let y = if y < 0.5 {
                    0.5
                } else if y > size as f64 + 0.5 {
                    size as f64 + 0.5
                } 
                else {
                    y
                };

                let i0 = x.floor() as usize;
                let i1 = i0 + 1;
                let j0 = y.floor() as usize;
                let j1 = j0 + 1;

                let s1 = x - i0 as f64;
                let s0 = 1.0 - s1;
                let t1 = y - j0 as f64;
                let t0 = 1.0 - t1;

                d[IX(i, j, size)] = s0 * (t0 * veloc_y[IX(i0, j0, size)] + t1 * veloc_y[IX(i0, j1, size)])
                                      + s1 * (t0 * veloc_y[IX(i1, j0, size)] + t1 * veloc_y[IX(i1, j1, size)]);
            }
        }
        Self::set_bnd(b, d, size);
    }
}

fn IX(x: usize, y: usize, width: usize) -> usize {
    x + y * width
}