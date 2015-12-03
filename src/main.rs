struct SiteIterator {
    cur: u64,
    max: (u64, u64, u64),
}


impl SiteIterator {

    pub fn new(shape: (u64, u64, u64)) -> SiteIterator {
        SiteIterator { cur: 0, max: shape }
    }

    fn cube(size: u64) ->  SiteIterator {
        SiteIterator {
            cur: 0,
            max: (size, size, size),
        }
    }
}


impl Iterator for SiteIterator {

    type Item = (u64, u64, u64);

    fn next(&mut self) ->  Option<(u64, u64, u64)> {
        if self.cur == self.max.0 * self.max.1 * self.max.2 {
            return None;
        }
        let x =  self.cur % self.max.0;
        let y = (self.cur / self.max.0) % self.max.1;
        let z = (self.cur / self.max.0 / self.max.1);
        self.cur += 1;
        Some((x, y, z))
    }
}

fn main() {
    println!("3x3x3 cube ->");
    for item in SiteIterator::cube(3) {
        println!("{:?}", item);
    }

    println!("1x3x3 layer ->");
    for item in SiteIterator::new((1, 3, 3)) {
        println!("{:?}", item);
    }
}
