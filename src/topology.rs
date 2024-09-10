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

    #[derive(Clone, PartialEq)]
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

    pub struct Sheaf {
        topology: Topology,
        sections: HashMap<OpenSet, Section>,
    }

    impl Sheaf {
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
            Sheaf { topology, sections }
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

        pub fn get_oset_from_section(&self, section: &Section) ->  Result<OpenSet, String> {
            if self.sections.iter().any(|(k, sec)| sec == section) == false {
                Err("Invalid section".to_string())
            } else {
            let mut oset: OpenSet = Vec::new();
            self.sections.iter().filter_map(|(k, _sec)| Some(oset = k.to_vec()));
            Ok(oset)
            }
        }

        pub fn restrict(&self, larger_oset: OpenSet, target_oset: &OpenSet) -> Result<Section, String> {
            if target_oset.iter().all(|point| larger_oset.contains(point)) == false {
                Err("Target Open Set is not a subset of the provided start set!".to_string())
            } else {
                let start = self.get_section(&larger_oset);
                let mut restricted_section: BTreeMap<LatticePoint, Observables> = BTreeMap::new();
                for (point, obs) in start.iter() {
                    restricted_section.insert(point.clone(), obs.clone());
                }
                Ok(restricted_section)
            }
        }

        pub fn can_glue(&self, first_section: &Section, second_section: &Section) -> bool {
            if self.sections.iter().any(|(k, sec)| sec == first_section || sec == second_section) == false  {
                false
            } else {
                let first_oset = self.get_oset_from_section(first_section).unwrap();
                let second_oset = self.get_oset_from_section(second_section).unwrap();
                if first_oset.iter().any(|point| second_oset.contains(point)) == false {
                    false
                } else {
                    let mut intersection: OpenSet = Vec::new();
                    let _ = first_oset.iter().filter(|point| 
                    if second_oset.iter().contains(point) {
                        intersection.push(point.to_vec());
                        Some(()).is_some()
                    } else {
                        None::<usize>.is_some()
                    }
                    );
                    let restricted_first = first_section.iter().filter_map(|(point, obs)| {
                        if intersection.iter().contains(point) {
                            Some((point.clone(), obs.clone()))
                        } else {
                            None
                        }
                    }).collect::<Section>();
                    let restricted_second = second_section.iter().filter_map(|(point, obs)| {
                        if intersection.iter().contains(point) {
                            Some((point.clone(), obs.clone()))
                        } else {
                            None
                        }
                    }).collect::<Section>();
                    if restricted_first.iter().zip(restricted_second).all(|((_, first), (_, second))| first == &second) == false {
                        false
                    } else {
                        true
                    }
                }
            }
        }

        pub fn glue(&self, first_section: &Section, second_section: &Section) -> Result<Section, String> {
            if self.can_glue(first_section, second_section) == false {
                Err("Cannot glue these sections".to_string())
            } else {
                
                Ok()
            }
        }            
    }


}
