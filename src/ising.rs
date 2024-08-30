use conv::prelude::*;
use itertools::Itertools;
use std::collections::HashMap;

use rand::Rng;

use crate::topology::*;

pub const BOLTZMANN: f64 = 1.380649e-23;

#[derive(Clone, Copy, PartialEq)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Clone)]
pub struct Lattice {
    pub dimension: usize,
    pub size: Vec<usize>,
}

impl Lattice {
    pub fn new(dimension: usize) -> Self {
        Lattice {
            dimension,
            size: Vec::new(),
        }
    }

    pub fn set_size(&mut self, size: Vec<usize>) {
        assert!(
            size.len() == self.dimension,
            "size vector does not match dimension of lattice"
        );
        self.size = size;
    }

    pub fn all_points(&self) -> impl Iterator<Item = LatticePoint> + '_ {
        (0..self.dimension)
            .map(|d| 0..self.size[d])
            .multi_cartesian_product()
    }
}

pub struct Ising {
    pub lattice: Lattice,
    pub spins: HashMap<Vec<usize>, Spin>,
    pub coupling: f64,
    pub applied_field: f64,
    pub temperature: f64,
    pub topology: Topology,
}

impl Ising {
    pub fn new(lattice: Lattice, coupling: f64, applied_field: f64, temperature: f64) -> Self {
        let spins = (0..lattice.dimension)
            .map(|d| 0..lattice.size[d])
            .multi_cartesian_product()
            .map(|idx| (idx, Spin::Up))
            .collect::<HashMap<Vec<usize>, Spin>>();
        let topology = Topology::new(lattice.clone());
        Ising {
            lattice,
            spins,
            coupling,
            applied_field,
            temperature,
            topology,
        }
    }

    pub fn get_spin(&self, idx: &[usize]) -> Result<Spin, &str> {
        if idx
            .iter()
            .zip(&self.lattice.size)
            .any(|(&i, &cap)| i >= cap)
        {
            return Err("Invalid Index");
        }
        Ok(self.spins.get(&idx.to_vec()).unwrap().clone())
    }

    pub fn set_spin(&mut self, idx: &[usize], spin: Spin) -> Result<(), &str> {
        if idx
            .iter()
            .zip(&self.lattice.size)
            .any(|(&i, &cap)| i >= cap)
        {
            return Err("Invalid Index");
        }
        self.spins.get(idx).replace(&spin);
        Ok(())
    }

    pub fn nearest_neighbor(&self, idx: &[usize]) -> Result<Vec<Vec<usize>>, &str> {
        if idx
            .iter()
            .zip(&self.lattice.size)
            .any(|(&i, &cap)| i >= cap)
        {
            return Err("Invalid Index");
        }
        let mut neighbors: Vec<Vec<usize>> = self
            .spins
            .keys()
            .filter(|&node| {
                node.iter()
                    .zip(idx)
                    .map(|(&n, &i)| abs_distance(n, i))
                    .sum::<usize>()
                    == 1
            })
            .cloned()
            .collect();
        Ok(neighbors)
    }

    pub fn local_energy(&self, idx: &[usize]) -> Result<f64, &str> {
        if idx
            .iter()
            .zip(&self.lattice.size)
            .any(|(&i, &cap)| i >= cap)
        {
            return Err("Invalid Index");
        }
        let local_spin = match self.get_spin(idx).unwrap() {
            Spin::Up => 1.0,
            Spin::Down => -1.0,
        };
        let field_energy = -self.applied_field * local_spin;
        let neighbor_energy: f64 = self
            .nearest_neighbor(idx)
            .unwrap()
            .iter()
            .map(|nidx| {
                let neighbor_spin = match self.spins.get(nidx).unwrap() {
                    Spin::Up => 1.0,
                    Spin::Down => -1.0,
                };
                -neighbor_spin * local_spin * self.coupling
            })
            .sum();
        Ok(field_energy + neighbor_energy)
    }

    pub fn total_energy(&self) -> f64 {
        self.spins
            .iter()
            .map(|(idx, _)| self.local_energy(idx).unwrap())
            .sum()
    }

    pub fn magnetization(&self) -> f64 {
        self.spins
            .iter()
            .map(|(_, &spin)| match spin {
                Spin::Up => 1.0,
                Spin::Down => -1.0,
            })
            .sum::<f64>()
            / self.spins.len().value_as::<f64>().unwrap()
    }

    pub fn metropolis_stepper(&mut self) {
        let mut rng = rand::thread_rng();
        let mut idx = Vec::new();
        for d in 0..self.lattice.dimension {
            idx.push(rng.gen_range(0..self.lattice.size[d]))
        }
        let init_energy = self.local_energy(idx.as_slice()).unwrap();
        let mut new_spin = match self.get_spin(idx.as_slice()).unwrap() {
            Spin::Up => Spin::Down,
            Spin::Down => Spin::Up,
        };
        self.set_spin(idx.as_slice(), new_spin);
        let energy_change = self.local_energy(idx.as_slice()).unwrap() - init_energy;
        if energy_change > 0.0
            && rng.gen::<f64>() > (-energy_change / (BOLTZMANN * self.temperature)).exp()
        {
        } else if energy_change < 0.0 {
        } else {
            new_spin = match new_spin {
                Spin::Up => Spin::Down,
                Spin::Down => Spin::Up,
            };
            self.set_spin(idx.as_slice(), new_spin);
        }
    }

    pub fn get_up_spin_set(&self) -> OpenSet {
        self.topology.open_set_from_spins(self, Spin::Up)
    }

    pub fn get_down_spin_set(&self) -> OpenSet {
        self.topology.open_set_from_spins(self, Spin::Down)
    }

    pub fn correlation(&self, idx: &[usize]) -> Result<f64, String> {
        if idx
            .iter()
            .zip(&self.lattice.size)
            .any(|(&i, &cap)| i >= cap)
        {
            return Err("Invalid Index".to_string());
        }
        let spin = match self.get_spin(idx).unwrap() {
            Spin::Up => 1.0,
            Spin::Down => -1.0,
        };
        let neighbors = self.nearest_neighbor(idx).unwrap();
        let neighbor_correlation = neighbors
            .iter()
            .map(|each| match self.get_spin(each.as_slice()).unwrap() {
                Spin::Up => 1.0 * spin,
                Spin::Down => -1.0 * spin,
            })
            .sum::<f64>()
            / neighbors.len().value_as::<f64>().unwrap();
        Ok(neighbor_correlation - self.magnetization().powf(2.0))
    }
}

pub fn abs_distance(a: usize, b: usize) -> usize {
    if a > b {
        a - b
    } else {
        b - a
    }
}
