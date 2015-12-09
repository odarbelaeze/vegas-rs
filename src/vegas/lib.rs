extern crate rand;

pub mod util;
pub mod lattice;
pub mod state;
pub mod energy;


#[macro_export]
macro_rules! hamiltonian {
    (@flatten $I: expr,) => (
        $I
        );
    (@flatten $I: expr, $J: expr, $($K:expr,)*) => (
        hamiltonian!(@flatten hamiltonian!($I, $J), $($K,)*)
        );
    ($I: expr) => (
        $I
        );
    ($I: expr, $J: expr) => (
        $crate::energy::ComposedEnergy::new($I, $J)
        );
    ($I: expr, $J: expr, $($K: expr),+) => (
        hamiltonian!(@flatten hamiltonian!($I, $J), $($K,)+)
        );
}
