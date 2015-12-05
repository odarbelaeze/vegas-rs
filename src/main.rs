use std::iter::Filter;
use std::slice::Iter;
use std::num;

mod sample;
use sample::util::super_mod;
use sample::lattice::Lattice;
use sample::lattice::Vertex;
use sample::lattice::Site;


struct Crystal {
    lims: Vec<usize>,
    nbhs: Vec<usize>,
}

impl Crystal {
    pub fn new(lattice: Lattice) -> Crystal
    {
        let mut lims = vec![0];
        let mut nbhs = vec![];
        for site in lattice.sites() {
            let pnbhs: Vec<usize>  = lattice.targets(&site).unwrap().iter().map(|site| {
                lattice.index(&site).unwrap()
            }).collect();
            let last = lims.last().unwrap().clone();
            lims.push(last + pnbhs.len());
            nbhs.extend(pnbhs.iter());
        }
        Crystal { lims: lims, nbhs: nbhs, }
    }

    fn nbhs_of<'a>(&'a self, item: usize) -> Option<&'a[usize]> {
        if item >= self.lims.len() - 1 {
            return None
        }
        let low = self.lims[item] as usize;
        let hi = self.lims[item + 1] as usize;
        Some(&self.nbhs[low..hi])
    }
}


fn main() {
    println!("test module {}", -11 % 10);
    assert_eq!(9, super_mod(-1,  10));
    assert_eq!(9, super_mod(-11, 10));
    assert_eq!(0, super_mod(-10, 10));

    let latt = Lattice {
        pbc: (true, true, false),
        shape: (10, 10, 10),
        natoms: 1,
        vertices: Vertex::list_for_cubic(),
    };

    match latt.inside(&Site { cell: (10, 10, 9), atom: 0 }) {
        None => panic!("La embarrasteis"),
        Some(site) => println!("Calidoso: {:?}", site),
    }

    match latt.inside(&Site { cell: (10, -1, 9), atom: 0 }) {
        None => panic!("La embarrasteis"),
        Some(site) => println!("Calidoso: {:?}", site),
    }

    match latt.inside(&Site { cell: (10, 10, 10), atom: 0 }) {
        None => println!("Calidoso!"),
        Some(site) => panic!("La embarrasteis {:?}", site),
    }

    match latt.inside(&Site { cell: (10, 10, 9), atom: 2 }) {
        None => println!("Calidoso!"),
        Some(site) => panic!("La embarrasteis {:?}", site),
    }

    println!("Neighbors attempt");
    let crystal = Crystal::new(latt);
    println!("len lims {}\nlen nbhs {}", crystal.lims.len(), crystal.nbhs.len());

    for item in &crystal.lims {
        println!("lims explicit {}", item);
    }

    for item in &crystal.nbhs {
        println!("lims explicit {}", item);
    }

    assert_eq!(vec![990, 909, 998, 989, 899], crystal.nbhs_of(999).unwrap());
}
