use linear_algebra::*;
use ndarray::arr2;
use num::rational::Ratio;

type T = i128;

fn main() {
  println!("Проверка задачи 1:");
  let a = arr2(&[[2, 1, 3, -1], [3, 2, 0, -2], [3, 1, 9, -1]])
    .map(|x: &i32| Ratio::<T>::from(T::from(x.clone())));
  let (r, c) = reduced_row_echelon_form(&a);
  println!("{}\n{}", r, c);
}
