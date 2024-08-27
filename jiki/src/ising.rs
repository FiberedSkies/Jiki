 use std::collections::HashMap;
 use itertools::Itertools;

 use rand::Rng;

pub const BOLTZMANN: f64 = 1.380649e-23;

pub enum Spin {
    Up,
    Down,
}

pub struct Lattice {
    dimension: usize,
    size: Vec<usize>,
}

impl Lattice {
    pub fn new(dimension: usize) -> Lattice {
        Lattice {
            dimension,
            size: Vec::new(),
        }
    }

    pub fn set_size(&mut self, size: Vec<usize>) {
        assert!(size.len() == self.dimension, "size vector does not match dimension of lattice");
        self.size = size;
    }
}

pub struct Ising {
    lattice: Lattice,
    spins: HashMap<Vec<usize>, Spin>,
    coupling: f64,
    applied_field: f64,
    temperature: f64,
}

impl Ising {
    pub fn new(lattice: Lattice, coupling: f64, applied_field: f64, temperature: f64) -> Ising {
        let spins = (0..lattice.dimension).map(|d| 0..lattice.size[d]).multi_cartesian_product().map(|idx| (idx, Spin::Up)).collect::<HashMap<Vec<usize>, Spin>>();
        Ising {
            lattice,
            spins,
            coupling,
            applied_field,
            temperature
        }
    }

    pub fn local_energy(&self, lattice_coords: Vec<usize>) -> Result<f64, &str> {
        
    }
}
