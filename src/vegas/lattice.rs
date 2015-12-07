//! Useful functions and data structures to build lattices

use util::super_mod;


#[derive(Debug)]
/// Represents a site, the cell represents where the site is
/// in the `Lattice` and the atom represents wich atom it is
/// within the unitcell
pub struct Site {
    cell: (i64, i64, i64),
    atom: u32,
}


/// Represents a lattice, it requires the periodicity of the lattice,
/// its shape, its number of atoms per site and the vertices as a list
/// of vertex descriptors
pub struct Lattice {
    pbc: (bool, bool, bool),
    shape: (u32, u32, u32),
    natoms: u32,
    vertices: Vec<Vertex>,
}


impl Lattice {

    /// Returns a site in the first image of a lattice according the lattice's
    /// periodicity, returns none if the site is outside the lattice
    pub fn inside(&self, site: &Site) ->  Option<Site> {

        if !(site.atom < self.natoms) {
            return None
        }

        let (mut x, mut y, mut z) = site.cell;
        let (sx, sy, sz) = (
            self.shape.0 as i64,
            self.shape.1 as i64,
            self.shape.2 as i64,
            );

        if !self.pbc.0  && (x < 0 || sx <= x) {
            return None
        } else {
            x = super_mod(x, sx);
        }

        if !self.pbc.1  && (y < 0 || sy <= y) {
            return None
        } else {
            y = super_mod(y, sy);
        }

        if !self.pbc.2  && (z < 0 || sz <= z) {
            return None
        } else {
            z = super_mod(z, sz);
        }

        Some(Site { cell: (x, y, z), atom: site.atom} )
    }

    pub fn sites(&self) -> SiteIterator {
        SiteIterator::new(self)
    }

    pub fn index(&self, site: &Site) -> Option<usize> {
        self.inside(site).map(|site| {
            let atom = site.atom as usize;
            let natm = self.natoms as usize;
            let (sx, sy) = (self.shape.0 as usize, self.shape.1 as usize);
            let (cx, cy, cz) = (
                site.cell.0 as usize,
                site.cell.1 as usize,
                site.cell.2 as usize,
                );
            natm * (sx * (sy * cz + cy) + cx) + atom
        })
    }

    pub fn targets(&self, site: &Site) -> Option<Vec<Site>> {
        self.inside(site).map(|site| {
            let mut tgts: Vec<Site> = vec![];
            for vx in &self.vertices {
                match vx.target_for(&site) {
                    None => continue,
                    Some(tgt) => {
                        match self.inside(&tgt) {
                            None => continue,
                            Some(tgt) => tgts.push(tgt),
                        }
                    },
                }
            }
            tgts
        })
    }

}


/// A little builder for lattices
pub struct LatticeBuilder {
    pbc: (bool, bool, bool),
    shape: (u32, u32, u32),
    natoms: u32,
    vertices: Vec<Vertex>,
}


/// A little builder for lattices
impl LatticeBuilder {
    pub fn new() -> LatticeBuilder {
        LatticeBuilder {
            pbc: (true, true, true),
            shape: (10u32, 10u32, 10u32),
            natoms: 1u32,
            vertices: Vec::new(),
        }
    }

    pub fn pbc(mut self, pbc: (bool, bool, bool)) -> LatticeBuilder {
        self.pbc = pbc;
        self
    }

    pub fn shape(mut self, shape: (u32, u32, u32)) -> LatticeBuilder {
        self.shape = shape;
        self
    }

    pub fn natoms(mut self, natoms: u32) -> LatticeBuilder {
        self.natoms = natoms;
        self
    }

    pub fn vertices(mut self, vertices: Vec<Vertex>) -> LatticeBuilder {
        self.vertices = vertices;
        self
    }

    pub fn finalize(self) -> Lattice {
        Lattice {
            pbc: self.pbc,
            shape: self.shape,
            natoms: self.natoms,
            vertices: self.vertices,
        }
    }
}


/// Iterates over the cells of a lattice
struct CellIterator {
    cur: u32,
    max: (u32, u32, u32),
}


impl CellIterator {
    pub fn new(lattice: &Lattice) -> CellIterator {
        CellIterator { cur: 0, max: lattice.shape }
    }
}


impl Iterator for CellIterator {

    type Item = (u32, u32, u32);

    fn next(&mut self) ->  Option<(u32, u32, u32)> {
        if self.cur == self.max.0 * self.max.1 * self.max.2 {
            return None;
        }
        let x =  self.cur % self.max.0;
        let y = (self.cur / self.max.0) % self.max.1;
        let z =  self.cur / self.max.0  / self.max.1 ;
        self.cur += 1;
        Some((x, y, z))
    }

}


/// Iterates over the sites of a lattice
pub struct SiteIterator {
    cell_it: CellIterator,
    cur_cell: Option<<CellIterator as Iterator>::Item>,
    cur_at: u32,
    max_at: u32,
}


impl SiteIterator {

    pub fn new(lattice: &Lattice) -> SiteIterator {
        let mut iter = CellIterator::new(lattice);
        SiteIterator {
            cur_cell: iter.next(),
            cell_it: iter,
            cur_at: 0,
            max_at: lattice.natoms,
        }
    }

}


impl Iterator for SiteIterator {

    type Item = Site;

    fn next(&mut self) ->  Option<Site> {
        if self.max_at == 0 {
            return None
        }
        if self.cur_at == self.max_at {
            self.cur_at = 0;
            self.cur_cell = self.cell_it.next();
        }
        let at = self.cur_at;
        self.cur_at = self.cur_at + 1;
        match self.cur_cell {
            None => None,
            Some((x, y, z)) => Some(Site {
                cell: (x as i64, y as i64, z as i64),
                atom: at,
            }),
        }
    }

}


/// Represents a vertex descriptor, for a vertex that can go beyond the
/// unit cell of a lattice.
pub struct Vertex {
    source: u32,
    target: u32,
    delta: (i64, i64, i64),
}


impl Vertex {

    fn target_for(&self, site: &Site) -> Option<Site> {
        if site.atom != self.source {
            return None
        }
        Some(Site {
            cell: (site.cell.0 + self.delta.0,
                   site.cell.1 + self.delta.1,
                   site.cell.2 + self.delta.2),
                   atom: self.target,
        })
    }

    pub fn list_for_cubic() ->  Vec<Vertex> {
        vec![
            Vertex{ source: 0, target: 0, delta: (1, 0, 0) },
            Vertex{ source: 0, target: 0, delta: (0, 1, 0) },
            Vertex{ source: 0, target: 0, delta: (0, 0, 1) },
            Vertex{ source: 0, target: 0, delta: (-1, 0, 0) },
            Vertex{ source: 0, target: 0, delta: (0, -1, 0) },
            Vertex{ source: 0, target: 0, delta: (0, 0, -1) },
        ]
    }

    pub fn list_for_hcp() ->  Vec<Vertex> {
        vec![
            // Zero in plane
            Vertex { source: 0, target: 0, delta: (1, 0, 0) },
            Vertex { source: 0, target: 0, delta: (0, 1, 0) },
            Vertex { source: 0, target: 0, delta: (1, 1, 0) },
            // Zero in plane backwards
            Vertex { source: 0, target: 0, delta: (-1,  0, 0) },
            Vertex { source: 0, target: 0, delta: ( 0, -1, 0) },
            Vertex { source: 0, target: 0, delta: (-1, -1, 0) },
            // One in plane
            Vertex { source: 1, target: 1, delta: (1, 0, 0) },
            Vertex { source: 1, target: 1, delta: (0, 1, 0) },
            Vertex { source: 1, target: 1, delta: (1, 1, 0) },
            // One in plane backwards
            Vertex { source: 1, target: 1, delta: (-1,  0, 0) },
            Vertex { source: 1, target: 1, delta: ( 0, -1, 0) },
            Vertex { source: 1, target: 1, delta: (-1, -1, 0) },
            // Zero with one
            Vertex { source: 0, target: 1, delta: ( 0,  0, 0) },
            Vertex { source: 0, target: 1, delta: (-1,  0, 0) },
            Vertex { source: 0, target: 1, delta: (-1, -1, 0) },
            // Zero with one downwards
            Vertex { source: 0, target: 1, delta: ( 0,  0, -1) },
            Vertex { source: 0, target: 1, delta: (-1,  0, -1) },
            Vertex { source: 0, target: 1, delta: (-1, -1, -1) },
            // One with zero
            Vertex { source: 1, target: 0, delta: ( 0,  0, 0) },
            Vertex { source: 1, target: 0, delta: ( 1,  0, 0) },
            Vertex { source: 1, target: 0, delta: ( 1,  1, 0) },
            // Zero with one upwards
            Vertex { source: 1, target: 0, delta: ( 0,  0,  1) },
            Vertex { source: 1, target: 0, delta: ( 1,  0,  1) },
            Vertex { source: 1, target: 0, delta: ( 1,  1,  1) },
        ]
    }

    pub fn list_for_honeycomb() ->  Vec<Vertex> {
        vec![
            Vertex{ source: 0, target: 1, delta: (0, 0, 0) },
            Vertex{ source: 1, target: 0, delta: (0, 0, 0) },
            Vertex{ source: 0, target: 1, delta: (1, 0, 0) },
            Vertex{ source: 1, target: 0, delta: (0, 1, 0) },
            Vertex{ source: 0, target: 1, delta: (-1, 0, 0) },
            Vertex{ source: 1, target: 0, delta: (0, -1, 0) },

            // For the 3D lulz
            Vertex{ source: 0, target: 0, delta: (0, 0, 1) },
            Vertex{ source: 1, target: 1, delta: (0, 0, 1) },
            Vertex{ source: 0, target: 0, delta: (0, 0, -1) },
            Vertex{ source: 1, target: 1, delta: (0, 0, -1) },
        ]
    }
}

// Tests

#[test]
fn testing_the_inside() {
    let latt = LatticeBuilder::new()
        .pbc((true, true, false))
        .vertices(Vertex::list_for_cubic())
        .finalize();
    assert!(latt.inside(&Site { cell: (10, 10, 9), atom: 0 }).is_some());
    assert!(latt.inside(&Site { cell: (10, -1, 9), atom: 0 }).is_some());
    assert!(latt.inside(&Site { cell: (10, 10, 10), atom: 0 }).is_none());
    assert!(latt.inside(&Site { cell: (10, 10, 9), atom: 2 }).is_none());
}
