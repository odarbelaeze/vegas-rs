//! Some utility functions that come, say, from another language
//! implementations.


/// Python-like module function, the negative values are re-interpreted as
/// positive values
///
/// # Examples
///
/// ```
/// assert_eq!(0, vegas::util::super_mod(-10, 10));
/// assert_eq!(9, vegas::util::super_mod(-1, 10));
/// ```
pub fn super_mod(num: i64, basis: i64) -> i64 {
    let num = num % basis;
    if num < 0 {
        return basis + num
    }
    num
}


#[test]
fn super_mod_works_as_in_python() {
    assert_eq!(0, super_mod(10, 10));
    assert_eq!(1, super_mod(1, 10));
    assert_eq!(1, super_mod(-9, 10));
    assert_eq!(9, super_mod(-1, 10));
}


#[test]
fn from_main() {
    assert_eq!(9, super_mod(-1,  10));
    assert_eq!(9, super_mod(-11, 10));
    assert_eq!(0, super_mod(-10, 10));
}
