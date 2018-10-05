extern crate nalgebra as na;
extern crate num;
use num::integer::{gcd, lcm};
use na::DMatrix;
type Mat = DMatrix<i64>;

#[derive(Debug, Clone)]
enum ElementaryOperation {
    Negate(usize),
    Switch(usize, usize),
    Add(usize, i64, usize),
}

fn elem_row_operate(m: &mut Mat, op: &ElementaryOperation) {
    match op {
        ElementaryOperation::Negate(row) => for x in m.row_mut(*row).iter_mut() {
            *x = -*x;
        },
        ElementaryOperation::Switch(row1, row2) => {
            m.swap_rows(*row1, *row2);
        }
        ElementaryOperation::Add(row1, c, row2) => {
            let cols = m.ncols();
            for i in 0..cols {
                m[(*row2, i)] = m[(*row2, i)] + *c * m[(*row1, i)];
            }
        }
    }
}

fn elem_col_operate(m: &mut Mat, op: &ElementaryOperation) {
    match op {
        ElementaryOperation::Negate(col) => for x in m.column_mut(*col).iter_mut() {
            *x = -*x;
        },
        ElementaryOperation::Switch(col1, col2) => {
            m.swap_columns(*col1, *col2);
        }
        ElementaryOperation::Add(col1, c, col2) => {
            let rows = m.nrows();
            for j in 0..rows {
                m[(j, *col2)] = m[(j, *col2)] + *c * m[(j, *col1)];
            }
        }
    }
}

fn row_elimination(m: &mut Mat, idx: usize) -> bool {
    let rows = m.nrows();
    let mut minval = 1000;
    let mut mini = 0;
    for i in idx..rows {
        let absval = m[(i, idx)].abs();
        if 0 < absval && absval < minval {
            minval = absval;
            mini = i;
        }
    }
    if minval == 1000 {
        return false;
    }
    elem_row_operate(m, &ElementaryOperation::Switch(idx, mini));
    if m[(idx, idx)] < 0 {
        elem_row_operate(m, &ElementaryOperation::Negate(idx));
    }
    let mut nz_idx = idx;
    for i in (idx + 1)..rows {
        let q = m[(i, idx)] / m[(idx, idx)];
        elem_row_operate(m, &ElementaryOperation::Add(idx, -q, i));
        if m[(i, idx)] != 0 {
            nz_idx = i;
        }
    }
    nz_idx != idx
}
fn col_elimination(m: &mut Mat, idx: usize) -> bool {
    let cols = m.ncols();
    let mut minval = 1000;
    let mut minj = 0;
    for j in idx..cols {
        let absval = m[(idx, j)].abs();
        if 0 < absval && absval < minval {
            minval = absval;
            minj = j;
        }
    }
    if minval == 1000 {
        return false;
    }
    elem_col_operate(m, &ElementaryOperation::Switch(idx, minj));
    if m[(idx, idx)] < 0 {
        elem_col_operate(m, &ElementaryOperation::Negate(idx));
    }
    let mut nz_idx = idx;
    for j in (idx + 1)..cols {
        let q = m[(idx, j)] / m[(idx, idx)];
        elem_col_operate(m, &ElementaryOperation::Add(idx, -q, j));
        if m[(idx, j)] != 0 {
            nz_idx = j;
        }
    }
    nz_idx != idx
}
fn diagonalize(m: &Mat) -> (Mat, usize) {
    let mut m = m.clone();
    let (rows, cols) = m.shape();
    let dsize = if rows < cols { rows } else { cols };
    for idx in 0..dsize {
        while row_elimination(&mut m, idx) {}
        while col_elimination(&mut m, idx) {}
    }
    let mut nzlen = 0;
    for idx in 0..dsize {
        if m[(idx, idx)].abs() > 0 {
            nzlen = nzlen + 1;
        }
    }
    (m, nzlen)
}

fn smith(m: &Mat) -> Vec<i64> {
    let (d, nzlen) = diagonalize(&m);
    let mut v = vec![0; nzlen];
    for i in 0..nzlen {
        v[i] = d[(i, i)];
    }
    convert_divisible_seq(&v)
}

fn is_divisible(v: &Vec<i64>) -> bool {
    let len = v.len();
    for i in 0..len - 1 {
        if v[i + 1] % v[i] != 0 {
            return false;
        }
    }
    true
}

fn convert_divisible_seq(v: &Vec<i64>) -> Vec<i64> {
    let mut v = v.clone();
    let len = v.len();
    while !is_divisible(&v) {
        v.sort();
        for i in 0..(len - 1) {
            for j in (i + 1)..len {
                let (g, l) = (gcd(v[i], v[j]), lcm(v[i], v[j]));
                v[i] = g;
                v[j] = l;
            }
        }
    }
    v
}

type Simplex = Vec<usize>;
type Chain = Vec<Vec<Simplex>>;

#[derive(Debug, Clone)]
struct SimplicialComplex {
    chain: Chain,
}

impl SimplicialComplex {
    fn from_chain(chain: &Chain) -> Self {
        SimplicialComplex {
            chain: chain.clone(),
        }
    }
    fn delta(simplex: &Simplex) -> Vec<(Simplex, i64)> {
        let mut result: Vec<(Simplex, i64)> = vec![(simplex.clone(), 0); simplex.len()];
        let mut sign = 1;
        for (idx, p) in &mut result.iter_mut().enumerate() {
            (*p).0.remove(idx);
            (*p).1 = sign;
            sign = -sign;
        }
        result
    }
    fn to_idx(&self, simplex: &Simplex) -> Option<usize> {
        let q = simplex.len() - 1;
        for (idx, s) in self.chain[q].iter().enumerate() {
            if s == simplex {
                return Some(idx);
            }
        }
        None
    }
    fn boundary_mat(&self, q: usize) -> Mat {
        let mut mat = Mat::from_element(self.chain[q].len(), self.chain[q + 1].len(), 0);
        for (j, s) in self.chain[q + 1].iter().enumerate() {
            let result: Vec<(Simplex, i64)> = SimplicialComplex::delta(&s);
            for (s, sign) in result {
                let i = self.to_idx(&s).expect("something wrong");
                mat[(i, j)] = sign;
            }
        }
        mat
    }
    fn homology_group(&self, q: usize) -> FinitelyGeneratedModule {
        let mat1 = self.boundary_mat(q);
        let mat2 = self.boundary_mat(q + 1);
        let (v1, v2) = (smith(&mat1), smith(&mat2));
        let rank = mat2.nrows() - v1.len() - v2.len();
        let mut torsion: Vec<i64> = Vec::new();
        for a in v2 {
            if a != 1 {
                torsion.push(a);
            }
        }
        FinitelyGeneratedModule {
            torsion: torsion,
            rank: rank,
        }
    }
}

struct FinitelyGeneratedModule {
    torsion: Vec<i64>,
    rank: usize,
}

impl FinitelyGeneratedModule {
    fn to_string(&self) -> String {
        let mut s = String::new();
        for a in &self.torsion {
            s = s + &format!("Z/{}Z⊕", a);
        }
        s + &format!("Z^{}", self.rank)
    }
}

fn main() {
    let c2 = vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]];
    let c1 = vec![
        vec![0, 1],
        vec![0, 2],
        vec![0, 3],
        vec![0, 4],
        vec![1, 2],
        vec![1, 3],
        vec![2, 3],
        vec![2, 4],
    ];
    let c0 = vec![vec![0], vec![1], vec![2], vec![3], vec![4]];
    let c = vec![c0, c1, c2];
    let sc = SimplicialComplex::from_chain(&c);
    println!("H_0 = {}", sc.homology_group(0).to_string());
}
