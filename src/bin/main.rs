extern crate vegas;

use vegas::lattice::Lattice;
use vegas::lattice::LatticeBuilder;
use vegas::lattice::Vertex;


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
            let pnbhs: Vec<usize>  = lattice.targets(&site)
                .expect("No sites ma frien").iter().map(|site| {
                    lattice.index(&site).expect("Site outside lattice")
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
    let latt = LatticeBuilder::new()
        .pbc((true, true, true))
        .shape((10, 10, 10))
        .vertices(Vertex::list_for_hcp())
        .natoms(2)
        .finalize();

    println!("Neighbors attempt");
    let crystal = Crystal::new(latt);
    println!("len lims {}\nlen nbhs {}", crystal.lims.len(), crystal.nbhs.len());
    // assert_eq!(vec![990, 909, 998, 989, 899], crystal.nbhs_of(999).unwrap());
}
