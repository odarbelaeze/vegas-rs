//! Most Hamiltonian models require pair interactions, a good way to express
//! these is through a graph using an adjacency matrix.

#[derive(Debug)]
struct SiteIterator {
    current: Vec<usize>,
    max: Vec<usize>,
    _carry: usize,
}

impl SiteIterator {
    fn new(max: Vec<usize>) -> SiteIterator {
        SiteIterator {
            current: vec![0; max.len()],
            max: max,
            _carry: 0,
        }
    }
}

impl Iterator for SiteIterator {

    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self._carry > 0 || self.max.iter().any(|m| *m == 0) || self.max.len() == 0 {
            None
        } else {
            let mut carry = 1;
            let old_current = self.current.clone();
            self.current = old_current
                .iter()
                .zip(self.max.iter())
                .map(|(c, m)| {
                    if c + carry >= *m {
                        carry = 1;
                        0
                    } else {
                        let tmp = c + carry;
                        carry = 0;
                        tmp
                    }
                }).collect();
            self._carry = carry;
            Some(old_current)
        }
    }
}

#[cfg(test)]
mod test{
    use super::SiteIterator;

    #[test]
    fn test_new() {
        let s = SiteIterator::new(vec![1, 1, 1]);
        assert!(s.current == vec![0, 0, 0]);
    }

    #[test]
    fn iter_over_zero_sized_yields_empty() {
        let s = SiteIterator::new(vec![0]);
        assert!(s.collect::<Vec<Vec<usize>>>().len() == 0);
    }

    #[test]
    fn simple_site_iterator() {
        let s = SiteIterator::new(vec![1]);
        assert!(s.collect::<Vec<Vec<usize>>>().len() == 1);
    }

    #[test]
    fn complex_site_iterator() {
        let s = SiteIterator::new(vec![10, 2]);
        let sites = s.collect::<Vec<Vec<usize>>>();
        assert!(sites.len() == 20);
    }

    #[test]
    fn more_complex_site_iterator() {
        let s = SiteIterator::new(vec![10, 10]);
        let sites = s.collect::<Vec<Vec<usize>>>();
        assert!(sites.len() == 100);
    }

    #[test]
    fn another_trivial_iterator() {
        let s = SiteIterator::new(vec![]);
        assert!(s.collect::<Vec<Vec<usize>>>().len() == 0);
    }
}
