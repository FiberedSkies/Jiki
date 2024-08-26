use std::collections::HashSet;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Clone, PartialEq, Debug)]
pub struct ExternalVars {
    temperature: f64,
    applied_field: (f64, f64),
}

#[derive(Clone)]
pub struct LocalState {
    env: ExternalVars,
    spin: Spin,
}

pub struct Ising2D {
    length: usize,
    width: usize,
    nodes: Vec<(usize, usize, LocalState)>,
    coupling: f64,
}

impl Ising2D {
    pub fn new(length: usize, width: usize, init_env: LocalState, coupling: f64) -> Self {
        let mut nodes: Vec<(usize, usize, LocalState)> = Vec::new();
        for i in 0..length {
            for j in 0..width {
                nodes.push((i, j, init_env.clone()));
            }
        }
        Self {
            length,
            width,
            nodes,
            coupling,
        }
    }
    pub fn get_state(&self, x: usize, y: usize) -> Result<&LocalState, &str> {
        if let Some(node) = self
            .nodes
            .iter()
            .filter(|&&(a, b, _)| a == x && b == y)
            .map(|(_, _, ls)| ls)
            .next()
        {
            Ok(node)
        } else {
            Err("Invalid Grid Coordinate")
        }
    }

    pub fn set_state(&mut self, x: usize, y: usize, new_state: LocalState) -> Result<(), &str> {
        if let Some(target_node) = self.nodes.iter().position(|&(a, b, _)| a == x && b == y) {
            let (a, b, c) = &mut self.nodes[target_node];
            *a = x;
            *b = y;
            *c = new_state;
            Ok(())
        } else {
            Err("Invalid Grid Coordinates")
        }
    }
    pub fn nearest_neighbor(&self, x: usize, y: usize) -> Vec<usize> {
        assert!(
            x < self.length && y < self.width,
            "Coordinates out of bounds!"
        );
        let mut neighbors: Vec<usize> = Vec::new();
        for node in &self.nodes {
            let ldist = abs_distance(node.0, x);
            let wdist = abs_distance(node.1, y);
            if (ldist == 1 && wdist == 0) || (ldist == 0 && wdist == 1) {
                let g = self
                    .nodes
                    .iter()
                    .position(|&(a, b, _)| a == node.0 && b == node.1)
                    .unwrap();
                neighbors.push(g);
            } else {
            }
        }
        neighbors
    }

    pub fn local_energy(&self, x: usize, y: usize) -> f64 {
        let state = self.get_state(x, y).unwrap();
        let spin = match state.spin {
            Spin::Up => 1.0,
            Spin::Down => -1.0,
        };
        let field_alignment = -state.env.applied_field.0 * state.env.applied_field.1 * spin;
        let neighbor_energy: f64 = self
            .nearest_neighbor(x, y)
            .iter()
            .map(|&idx| {
                let neighbor_state = &self.nodes[idx].2;
                let neighbor_spin = match neighbor_state.spin {
                    Spin::Up => 1.0,
                    Spin::Down => -1.0,
                };
                -self.coupling * spin * neighbor_spin
            })
            .sum();
        field_alignment + neighbor_energy
    }

    pub fn total_energy(&self) -> f64 {
        (0..self.length)
            .flat_map(|x| (0..self.width).map(move |y| self.local_energy(x, y)))
            .sum::<f64>()
            / 2.0
    }

    pub fn net_magnetization(&self) -> f64 {
        (0..self.length)
            .flat_map(|x| {
                (0..self.width).map(move |y| match self.get_state(x, y).unwrap().spin {
                    Spin::Up => 1.0,
                    Spin::Down => -1.0,
                })
            })
            .sum::<f64>()
            / 2.0
    }
}

pub fn abs_distance(a: usize, b: usize) -> usize {
    if a > b {
        a - b
    } else {
        b - a
    }
}
