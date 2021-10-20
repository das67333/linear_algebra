use linear_algebra::*;

fn timer(m: Matrix) -> rug::Integer {
    let calculation_timer = std::time::Instant::now();
    let det = m.determinant_bareiss();
    println!(
        "{:?} {} ms",
        m.size(),
        calculation_timer.elapsed().as_millis()
    );
    det
}

fn main() {
    let m;
    {
        let loading_timer = std::time::Instant::now();
        m = Matrix::from_file("matrix.csv").expect("loading matrix failed");
        println!("{} ms", loading_timer.elapsed().as_millis());
    }
    println!("{}", m);
    for n in 2..10 {
        println!("{}", timer(m.reduced(n * 50, n * 50)).to_string().len());
    }
}

/*
implemented:
    matr.size()
    Matrix::new(vec![vec![a, b], vec![c, d]])
    Matrix::from_file("matrix.csv")
    Matrix::zero(n, m)
    Matrix::id(n)
    matr.transposed()
    matrix![[a, b], [c, d]]
    matr[i][j]
    matr.to_string()
    matr *= num
    matr * num
    num * matr
    matr += matr
    matr + matr
    matr *= matr
    matr * matr
*/
