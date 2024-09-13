use itertools::Itertools;

use crate::ising::*;

use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

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
        (0..lattice.dimension)
            .map(|d| 0..lattice.size[d])
            .multi_cartesian_product()
            .collect::<Vec<Vec<usize>>>().iter().map(|p| {
                basis.insert(vec![p.to_vec()])
            });
        Topology { lattice, basis }
    }

    pub fn intersection(&self, mut sets: Vec<OpenSet>) -> OpenSet {
        if sets.is_empty() {
            return Vec::new()
        }
        let mut intersection = sets.pop().unwrap();
        for set in sets {
            intersection = intersection.into_iter().filter(|&point| set.contains(&point)).collect();
        };
        intersection
    }

    pub fn union(&self, sets: Vec<OpenSet>) -> OpenSet {
        if sets.is_empty() {
            return Vec::new()
        }
        let mut result = Vec::new();
        for set in sets {
            result.extend(set.iter().cloned());
        }
        result.into_iter().collect()
    }

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
    pub enum Observable {
        Energy,
        Spin,
        Correlation,
    }

    impl Observable {
        pub fn compute(ising: &Ising, idx: LatticePoint, obs: Observable) -> Result<f64, String> {
            if idx
                .iter()
                .zip(&ising.lattice.size)
                .any(|(&i, &cap)| i >= cap)
            {
                return Err("Invalid Index".to_string());
            }
            let result = match obs {
                Observable::Energy => ising.local_energy(idx.as_slice()).unwrap(),
                Observable::Spin => match ising.get_spin(idx.as_slice()).unwrap() {
                    Spin::Up => 1.0,
                    Spin::Down => -1.0,
                },
                Observable::Correlation => ising.correlation(idx.as_slice()).unwrap(),
            };
            Ok(result)
        }
    }

    type Section = BTreeMap<LatticePoint, f64>;

    pub struct Sheaf {
        topology: Topology,
        energy_sections: HashMap<OpenSet, Section>,
        spin_sections: HashMap<OpenSet, Section>,
        correlation_sections: HashMap<OpenSet, Section>,
    }

    impl Sheaf {
        pub fn new(topology: Topology, ising: &Ising) -> Self {
            let mut energy_sections = HashMap::new();
            let mut spin_sections = HashMap::new();
            let mut correlation_sections = HashMap::new();
            topology.basis.iter().for_each(|openset| {
                let mut energy_section = BTreeMap::new();
                let mut spin_section = BTreeMap::new();
                let mut correlation_section = BTreeMap::new();
                openset.iter().for_each(|point| {
                    energy_section.insert(
                        point.clone(),
                        Observable::compute(ising, point.clone(), Observable::Energy).unwrap(),
                    );
                    spin_section.insert(
                        point.clone(),
                        Observable::compute(ising, point.clone(), Observable::Spin).unwrap(),
                    );
                    correlation_section.insert(
                        point.clone(),
                        Observable::compute(ising, point.clone(), Observable::Correlation).unwrap(),
                    );
                });
                energy_sections.insert(openset.clone(), energy_section);
                spin_sections.insert(openset.clone(), spin_section);
                correlation_sections.insert(openset.clone(), correlation_section);
            });
            Sheaf { topology, energy_sections, spin_sections, correlation_sections }
        }

        pub fn get_section(&self, open_set: &OpenSet, obs: Observable) -> Section {
            let mut section = BTreeMap::new();
            for point in open_set {
                match obs {
                    Observable::Energy => {if let Some((_, basis_section)) = self
                        .energy_sections
                        .iter()
                        .find(|(basis_set, _)| basis_set.contains(point))
                    {
                        if let Some(observables) = basis_section.get(point) {
                            section.insert(point.clone(), observables.clone());
                        }
                    }},
                    Observable::Spin => {if let Some((_, basis_section)) = self
                        .spin_sections
                        .iter()
                        .find(|(basis_set, _)| basis_set.contains(point))
                    {
                        if let Some(observables) = basis_section.get(point) {
                            section.insert(point.clone(), observables.clone());
                        }
                    }},
                    Observable::Correlation => {if let Some((_, basis_section)) = self
                        .correlation_sections
                        .iter()
                        .find(|(basis_set, _)| basis_set.contains(point))
                    {
                        if let Some(observables) = basis_section.get(point) {
                            section.insert(point.clone(), observables.clone());
                        }
                    }},
                }
            }
            section
        }

        pub fn get_oset_from_section(&self, section: &Section, obs: Observable) ->  Result<OpenSet, String> {
            match obs {
                Observable::Energy => {if self.energy_sections.iter().any(|(k, sec)| sec == section) == false {
                    Err("Invalid section".to_string())
                } else {
                let mut oset: OpenSet = Vec::new();
                self.energy_sections.iter().filter_map(|(k, _sec)| Some(oset = k.to_vec()));
                Ok(oset)
                }},
                Observable::Spin => {if self.spin_sections.iter().any(|(k, sec)| sec == section) == false {
                    Err("Invalid section".to_string())
                } else {
                let mut oset: OpenSet = Vec::new();
                self.spin_sections.iter().filter_map(|(k, _sec)| Some(oset = k.to_vec()));
                Ok(oset)
                }},
                Observable::Correlation => {if self.correlation_sections.iter().any(|(k, sec)| sec == section) == false {
                    Err("Invalid section".to_string())
                } else {
                let mut oset: OpenSet = Vec::new();
                self.correlation_sections.iter().filter_map(|(k, _sec)| Some(oset = k.to_vec()));
                Ok(oset)
                }},
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

        pub fn glue(&self, first_section: &Section, second_section: &Section) -> Result<Section, String> {
            if self.sections.iter().any(|(k, sec)| sec == first_section || sec == second_section) == false  {
                Err("Invalid sections!".to_string())
            } else {
                let first_oset = self.get_oset_from_section(first_section).unwrap();
                let second_oset = self.get_oset_from_section(second_section).unwrap();
                if first_oset.iter().any(|point| second_oset.contains(point)) == false {
                    Err("The open sets do not overlap".to_string())
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
                        Err("the sections do not agree on the overlap".to_string())
                    } else {
                        let mut final_section: Section = BTreeMap::new();
                        let remaining_second = second_section.iter().filter_map(|(point, obs)| {
                            if intersection.iter().contains(point) {
                                None
                            } else {
                                Some((point.clone(), obs.clone()))
                            }
                        }).collect::<Section>();
                        let remaining_first = first_section.iter().filter_map(|(point, obs)| {
                            if intersection.iter().contains(point) {
                                None
                            } else {
                                Some((point.clone(), obs.clone()))
                            }
                        }).collect::<Section>();
                        final_section.extend(remaining_first);
                        final_section.extend(remaining_second);
                        final_section.extend(restricted_first);
                        Ok(final_section)
                    }
                }
            }
        }            
    }


}
