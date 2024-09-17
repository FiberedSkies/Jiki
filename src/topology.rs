use itertools::Itertools;

use crate::ising::*;

use std::collections::HashSet;

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
        let _ = (0..lattice.dimension)
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
            intersection = intersection.into_iter().filter(|point| set.contains(point)).collect();
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

    #[derive(Clone, PartialEq, Eq, Hash)]
    pub enum Observable {
        Energy,
        Spin,
        Correlation,
    }

    impl Observable {
        pub fn compute(ising: &Ising, idx: &LatticePoint, obs: Observable) -> Result<f64, String> {
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

    type Section<'a> = BTreeMap<&'a LatticePoint, f64>;

    pub struct Sheaf<'a> {
        topology: &'a Topology,
        sections: HashMap<&'a Observable, HashMap<&'a OpenSet, Section<'a>>>
    }

    impl<'a> Sheaf<'a> {
        pub fn new(topology: &'a Topology, ising: &Ising) -> Self {
            let mut all_sections = HashMap::new();
            for obs in &[Observable::Energy, Observable::Spin, Observable::Correlation] {
                let mut obs_sections = HashMap::new();
                for oset in &topology.basis {
                    let section: Section = oset.iter().map(|point| {
                        (point, Observable::compute(ising, point, obs.clone()).unwrap())
                    }).collect();
                    obs_sections.insert(oset, section);
                }
                all_sections.insert(obs, obs_sections);
            }
            Sheaf { topology , sections: all_sections }
        }

        pub fn get_sections(&mut self, open_set:&'a OpenSet) -> Vec<&Section<'a>> {
            let mut secs = Vec::new();
            for obs in &[Observable::Energy, Observable::Spin, Observable::Correlation] {
                let mut obs_section_over_oset: Section = BTreeMap::new();
                for point in open_set {
                    if let Some((_, sections)) = self.sections.get(obs).unwrap().iter().find(|(basis, _)|basis.contains(&point)) {
                        obs_section_over_oset.insert(&point, sections.get(&point).unwrap().clone());
                    }
                }
                self.sections.get_mut(obs).unwrap().insert(&open_set, obs_section_over_oset);
            }
            for obs in &[Observable::Energy, Observable::Spin, Observable::Correlation] {
                secs.push(self.sections.get(obs).unwrap().get(open_set).unwrap());
            }
            secs
        }

        pub fn restrict_sections(&mut self, open_set:&'a OpenSet, smaller_set: &'a OpenSet) -> Result<Vec<Section<'a>>, String> {
            if smaller_set.iter().all(|point| open_set.contains(point)) == false {
                Err("Target Open Set is not a subset of the provided start set!".to_string())
            } else {
                let initial_sections = self.get_sections(open_set);
                let mut restricted_sections = Vec::<Section<'a>>::new();
                for sec in initial_sections {
                    let mut restricted_sec = BTreeMap::new();
                    for point in smaller_set {
                        let val  = sec.iter().find_map(|(&point, obs)| {
                            if smaller_set.contains(point) {
                                Some(obs.clone())
                            } else {
                                None
                            }
                        }).unwrap();
                        restricted_sec.insert(point, val);
                    }
                    restricted_sections.push(restricted_sec);
                }
                for (obs, section) in [Observable::Energy, Observable::Spin, Observable::Correlation].iter().zip(restricted_sections.clone()) {
                    self.sections.get_mut(obs).unwrap().insert(smaller_set, section);
                }
                Ok(restricted_sections)
            }
        }

        pub fn glue(&mut self, open_sets: &'a mut Vec<OpenSet>) -> Result<Vec<Section<'a>>, String> {
            let intersection = self.topology.intersection(open_sets.clone());
            if intersection.len() == 0 {
                return Err("Open sets provided do not overlap!".to_string())
            }
            let mut all_sections: Vec<Vec<&Section>> = Vec::new();
            for oset in open_sets {
                all_sections.push(self.sections.values().map(|obs_sections| {
                    obs_sections.get(oset).unwrap()
                }).collect());
            }
            let mut all_restricted_sections: Vec<Vec<Section>> = all_sections.iter().map(|sections| {
                let restricted_secs: Vec<Section> = sections.iter().map(|&sections| {
                    let rest_sec: Section = sections.iter().filter_map(|(&point, &val)| {
                        if intersection.contains(point) {
                            Some((point, val))
                        } else {
                            None
                        }
                    }).collect();
                    rest_sec
                }).collect();
                restricted_secs
            }).collect();

            let mut compliment_sections: Vec<Vec<Section>> = all_sections.iter().map(|sections| {
                let restricted_secs: Vec<Section> = sections.iter().map(|&sections| {
                    let rest_sec: Section = sections.iter().filter_map(|(&point, &val)| {
                        if !intersection.contains(point) {
                            Some((point, val))
                        } else {
                            None
                        }
                    }).collect();
                    rest_sec
                }).collect();
                restricted_secs
            }).collect();

            let mut glued_sections: Vec<Section> = Vec::new();
            let mut comparison = all_restricted_sections.pop().unwrap();
            for obs in 0..3 {
                let mut glued_observable: Section = BTreeMap::new(); 
                let mut can_glue = Vec::new();
                for point in &intersection {
                    let comparison_val = comparison[obs].get(point);
                    can_glue.push(all_restricted_sections.iter().all(|sections| sections[obs].get(point) == comparison_val));
                }
                if can_glue.into_iter().all(|val| val == true) {
                    glued_observable.append(&mut comparison[obs]);
                    while let Some(sec) = compliment_sections.pop() {
                        let mut obs_sec = sec[obs].clone();
                        glued_observable.append(&mut obs_sec);
                    }
                } else {}
                glued_sections.push(glued_observable);
            }
            Ok(glued_sections)
        }         
    }


}
