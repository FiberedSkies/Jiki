use itertools::Itertools;

use crate::ising::*;

use std::collections::{BTreeSet, HashMap, HashSet};

pub type LatticePoint = Vec<usize>;
pub type OpenSet = Vec<LatticePoint>;

pub struct Topology {
    lattice: Lattice,
    basis: HashSet<OpenSet>,
}

impl Topology {
    pub fn new(lattice: Lattice) -> Self {
        let mut basis = HashSet::new();
        basis.insert(Vec::new());
        basis.insert(
            (0..lattice.dimension)
                .map(|d| 0..lattice.size[d])
                .multi_cartesian_product()
                .collect::<Vec<Vec<usize>>>(),
        );
        Topology { lattice, basis }
    }

    pub fn add_basis(&mut self, set: OpenSet) {
        self.basis.insert(set);
    }

    pub fn union<'a>(&self, sets: impl IntoIterator<Item = &'a OpenSet>) -> OpenSet {
        let mut result = BTreeSet::new();
        for set in sets {
            result.extend(set.iter().cloned());
        }
        result.into_iter().collect()
    }

    pub fn is_open(&self, set: &OpenSet) -> bool {
        let set: BTreeSet<_> = set.iter().cloned().collect();
        self.basis.iter().all(|base| {
            base.iter().all(|point| set.contains(point))
                || base.iter().all(|point| !set.contains(point))
        })
    }

    pub fn is_closed(&self, set: &OpenSet) -> bool {
        let complement: OpenSet = self
            .lattice
            .all_points()
            .filter(|point| !set.contains(point))
            .collect();
        self.is_open(&complement)
    }

    pub fn get_open_neighborhood(&self, point: &LatticePoint) -> OpenSet {
        self.basis
            .iter()
            .filter(|set| set.contains(point))
            .fold(Vec::new(), |acc, set| {
                self.union([&acc, set].iter().copied())
            })
    }

    pub fn closure(&self, set: &OpenSet) -> OpenSet {
        self.lattice
            .all_points()
            .filter(|point| {
                let neighborhood = self.get_open_neighborhood(point);
                !neighborhood.iter().all(|p| !set.contains(p))
            })
            .collect()
    }
}

impl Topology {
    pub fn open_set_from_spins(&self, ising: &Ising, spin: Spin) -> OpenSet {
        ising
            .lattice
            .all_points()
            .filter(|point| ising.get_spin(point).unwrap() == spin)
            .collect()
    }
}

pub mod sheaf {
    use std::collections::{BTreeMap, HashMap};

    use super::*;

    #[derive(Clone)]
    pub struct Observables {
        energy: f64,
        spin: f64,
        correlation: f64,
    }

    impl Observables {
        pub fn compute(ising: &Ising, idx: LatticePoint) -> Result<Self, String> {
            if idx
                .iter()
                .zip(&ising.lattice.size)
                .any(|(&i, &cap)| i >= cap)
            {
                return Err("Invalid Index".to_string());
            }
            Ok(Observables {
                energy: ising.local_energy(idx.as_slice()).unwrap(),
                spin: match ising.get_spin(idx.as_slice()).unwrap() {
                    Spin::Up => 1.0,
                    Spin::Down => -1.0,
                },
                correlation: ising.correlation(idx.as_slice()).unwrap(),
            })
        }
    }

    type Section = BTreeMap<LatticePoint, Observables>;

    pub struct PreSheaf {
        topology: Topology,
        sections: HashMap<OpenSet, Section>,
    }

    impl PreSheaf {
        pub fn new(topology: Topology, ising: &Ising) -> Self {
            let mut sections = HashMap::new();
            topology.basis.iter().for_each(|openset| {
                let mut section = BTreeMap::new();
                openset.iter().for_each(|point| {
                    section.insert(
                        point.clone(),
                        Observables::compute(ising, point.clone()).unwrap(),
                    );
                });
                sections.insert(openset.clone(), section);
            });
            PreSheaf { topology, sections }
        }

        pub fn get_section(&self, open_set: &OpenSet) -> Section {
            let mut section = BTreeMap::new();
            for point in open_set {
                if let Some((_, basis_section)) = self
                    .sections
                    .iter()
                    .find(|(basis_set, _)| basis_set.contains(point))
                {
                    if let Some(observables) = basis_section.get(point) {
                        section.insert(point.clone(), observables.clone());
                    }
                }
            }
            section
        }
    }
}
