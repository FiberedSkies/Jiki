use std::collections::HashSet;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Clone, PartialEq, Debug)]
pub struct ExternalVars {
    temperature: f64,
    applied_field: Vec<f64>,
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
                nodes.push((i,j,init_env.clone()));
            }
        }
        Self {
            length,
            width,
            nodes,
            coupling,
        }
    }
    pub fn get_state(&self, x: usize, y: usize) -> Option<&LocalState> {
        let node = self.nodes.iter().filter(|&&(a, b, _)| a == x && b == y).map(|(_,_,ls)| ls).next();
        node
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
}