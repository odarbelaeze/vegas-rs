//! Some utility functions that come, say, from another language
//! implementations.


/// Python-like module function, the negative values are re-interpreted as
/// positive values
///
/// # Examples
///
/// ```
/// use main::super_mod;
/// assert_eq!(0, super_mod(-10, 10));
/// ```
pub fn super_mod(num: i64, basis: i64) ->  i64 {
    let num = num % basis;
    if num < 0 {
        return basis + num
    }
    num
}
