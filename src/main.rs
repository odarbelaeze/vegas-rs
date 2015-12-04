use std::iter::Filter;
use std::slice::Iter;
use std::num;


#[derive(Debug)]
struct Site {
    cell: (i64, i64, i64),
    atom: u32,
}


struct CellIterator {
    cur: i64,
    max: (i64, i64, i64),
}


struct Lattice {
    pbc: (bool, bool, bool),
    shape: (i64, i64, i64),
    natoms: u32,
    neighbors: NeighborList,
}


fn super_mod(num: i64, basis: i64) ->  i64 {
    let num = num % basis;
    if num < 0 {
        return basis + num
    }
    num
}


impl Lattice {
    fn inside(&self, site: &Site) ->  Option<Site> {
        if !(site.atom < self.natoms) {
            return None
        }

        let (mut x, mut y, mut z) = site.cell;

        if (!self.pbc.0  && (x < 0 || self.shape.0 <= x)) {
            return None
        } else {
            x = super_mod(x, self.shape.0);
        }

        if (!self.pbc.1  && (y < 0 || self.shape.1 <= y)) {
            return None
        } else {
            y = super_mod(y, self.shape.1);
        }

        if (!self.pbc.2  && (z < 0 || self.shape.2 <= z)) {
            return None
        } else {
            z = super_mod(z, self.shape.2);
        }

        Some(Site { cell: (x, y, z), atom: site.atom} )
    }
}


enum Template {
    cubic, hcp,
}


impl CellIterator {

    pub fn new(shape: (i64, i64, i64)) -> CellIterator {
        CellIterator { cur: 0, max: shape }
    }

    fn cube(size: i64) ->  CellIterator {
        CellIterator {
            cur: 0,
            max: (size, size, size),
        }
    }

    fn shape(&self) ->  (i64, i64, i64) {
        self.max
    }
}


impl Iterator for CellIterator {

    type Item = (i64, i64, i64);

    fn next(&mut self) ->  Option<(i64, i64, i64)> {
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


struct BasisIterator {
    cell_it: CellIterator,
    cur_cell: Option<<CellIterator as Iterator>::Item>,
    cur_at: u32,
    max_at: u32,
}


impl BasisIterator {

    pub fn new(iter: CellIterator) -> BasisIterator {
        let mut iter = iter;
        BasisIterator {
            cur_cell: iter.next(),
            cell_it: iter,
            cur_at: 0,
            max_at: 1,
        }
    }

    pub fn new_with_atoms(iter: CellIterator, natoms: u32) -> BasisIterator {
        let mut iter = iter;
        BasisIterator {
            cur_cell: iter.next(),
            cell_it: iter,
            cur_at: 0,
            max_at: natoms,
        }
    }

    fn natoms(&self) -> u32 {
        self.max_at
    }
}


impl Iterator for BasisIterator {

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
            Some(cell) => Some(Site { cell: cell, atom: at }),
        }
    }
}


struct Neighbor {
    source: u32,
    target: u32,
    delta: (i64, i64, i64),
}


impl Neighbor {
    pub fn new(source: u32, target: u32, delta: (i64, i64, i64)) -> Neighbor {
        Neighbor { source: source, target: target, delta: delta }
    }

    pub fn of(&self, site: &Site) ->  Option<Site> {
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
}


struct NeighborList {
    list: Vec<Neighbor>,
}

impl NeighborList {

    fn cubic() ->  NeighborList {
        NeighborList {
            list: vec![
                Neighbor::new(0, 0, (1, 0, 0)),
                Neighbor::new(0, 0, (0, 1, 0)),
                Neighbor::new(0, 0, (0, 0, 1)),
                Neighbor::new(0, 0, (-1, 0, 0)),
                Neighbor::new(0, 0, (0, -1, 0)),
                Neighbor::new(0, 0, (0, 0, -1)),
            ],
        }
    }

    fn hcp() ->  NeighborList {
        NeighborList {
            list: vec![
                // Zero in plane
                Neighbor::new(0, 0, (1, 0, 0)),
                Neighbor::new(0, 0, (0, 1, 0)),
                Neighbor::new(0, 0, (1, 1, 0)),
                // Zero in plane backwards
                Neighbor::new(0, 0, (-1, 0, 0)),
                Neighbor::new(0, 0, (0, -1, 0)),
                Neighbor::new(0, 0, (-1, -1, 0)),
                // One in plane
                Neighbor::new(1, 1, (1, 0, 0)),
                Neighbor::new(1, 1, (0, 1, 0)),
                Neighbor::new(1, 1, (1, 1, 0)),
                // One in plane backwards
                Neighbor::new(1, 1, (-1, 0, 0)),
                Neighbor::new(1, 1, (0, -1, 0)),
                Neighbor::new(1, 1, (-1, -1, 0)),
                // Zero with one
                Neighbor::new(0, 1, ( 0,  0, 0)),
                Neighbor::new(0, 1, (-1,  0, 0)),
                Neighbor::new(0, 1, (-1, -1, 0)),
                // Zero with one downwards
                Neighbor::new(0, 1, ( 0,  0, -1)),
                Neighbor::new(0, 1, (-1,  0, -1)),
                Neighbor::new(0, 1, (-1, -1, -1)),
                // One with zero
                Neighbor::new(1, 0, ( 0,  0, 0)),
                Neighbor::new(1, 0, ( 1,  0, 0)),
                Neighbor::new(1, 0, ( 1,  1, 0)),
                // Zero with one upwards
                Neighbor::new(1, 0, ( 0,  0,  1)),
                Neighbor::new(1, 0, ( 1,  0,  1)),
                Neighbor::new(1, 0, ( 1,  1,  1)),
            ],
        }
    }
}


struct Crystal {
    lims: Vec<usize>,
    nbhs: Vec<usize>,
}

impl Crystal {
    pub fn new(
        neighbors: NeighborList,
        basis: BasisIterator,
        pbc: (bool, bool, bool)
        ) -> Crystal
    {
        let mut lims = vec![0];
        let mut nbhs = vec![];
        let shape = basis.cell_it.shape();
        let natoms = basis.natoms();
        for site in basis {
            let mut pnbhs = vec![];
            for nbh in neighbors.list.iter() {
                match nbh.of(&site) {
                    Some(site) => { pnbhs.push(0); }
                    None => {}
                }
            }
            let last = lims.last().unwrap().clone();
            lims.push(last + pnbhs.len());
            nbhs.extend(pnbhs.iter());
        }
        Crystal { lims: lims, nbhs: nbhs, }
    }
}


fn main() {
    println!("3x3x3 cube ->");
    for item in CellIterator::cube(3) {
        println!("{:?}", item);
    }

    println!("1x3x3 layer ->");
    for item in CellIterator::new((1, 3, 3)) {
        println!("{:?}", item);
    }

    println!("Cell with atoms");
    for site in BasisIterator::new(CellIterator::cube(2)) {
        println!("{:?}", site);
    }

    println!("Cell with multiple atoms");
    for site in BasisIterator::new_with_atoms(CellIterator::cube(2), 2) {
        println!("{:?}", site);
    }

    println!("Cell with multiple atoms");
    let upper = Neighbor::new(0, 0, (0, 0, 1));
    for site in BasisIterator::new_with_atoms(CellIterator::cube(2), 2) {
        let sup = upper.of(&site);
        println!("{:?} -> upper (of type 0) -> {:?}", site, sup);
    }

    println!("Neighbors attempt");
    let crystal = Crystal::new(
        NeighborList::hcp(),
        BasisIterator::new_with_atoms(CellIterator::cube(3), 2),
        (true, true, true));
    println!("len lims {}\nlen nbhs {}", crystal.lims.len(), crystal.nbhs.len());
    for item in crystal.lims {
        println!("lims explicit {}", item);
    }

    println!("test module {}", -11 % 10);
    assert_eq!(9, super_mod(-1,  10));
    assert_eq!(9, super_mod(-11, 10));
    assert_eq!(0, super_mod(-10, 10));

    let latt = Lattice {
        pbc: (true, true, false),
        shape: (10, 10, 10),
        natoms: 2,
        neighbors: NeighborList::cubic(),
    };

    match latt.inside(&Site { cell: (10, 10, 9), atom: 1 }) {
        None => panic!("La embarrasteis"),
        Some(site) => println!("Calidoso: {:?}", site),
    }
}
