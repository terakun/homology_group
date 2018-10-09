extern crate nalgebra as na;
extern crate num;
use num::integer::{gcd, lcm};
use na::DMatrix;
use std::collections::BTreeSet;
use std::io::BufReader;
use std::io::prelude::*;
use std::fs::File;
use std::env;
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
    // 0でない絶対値最小の要素の添字を求める
    let min_i = match m.column(idx)
        .iter()
        .enumerate()
        .skip(idx)
        .map(|(i, x)| (x.abs(), i))
        .filter(|(absx, _)| *absx > 0)
        .min()
    {
        Some((_, i)) => i,
        None => {
            return false;
        }
    };
    elem_row_operate(m, &ElementaryOperation::Switch(idx, min_i));
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
    // 0でない絶対値最小の要素の添字を求める
    let min_j = match m.row(idx)
        .iter()
        .enumerate()
        .skip(idx)
        .map(|(j, x)| (x.abs(), j))
        .filter(|(absx, _)| *absx > 0)
        .min()
    {
        Some((_, j)) => j,
        None => {
            return false;
        }
    };
    elem_col_operate(m, &ElementaryOperation::Switch(idx, min_j));
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
    let mut v = Vec::new();
    let mut i = 0;
    while v.len() < nzlen {
        if d[(i, i)] > 0 {
            v.push(d[(i, i)]);
        }
        i = i + 1;
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
    fn from_simplices(simplices: &Vec<Simplex>) -> Self {
        let dim = simplices[0].len() - 1;
        let mut chain: Chain = vec![Vec::new(); dim + 1];
        let mut simplices = simplices.clone();
        for simplex in simplices.iter_mut() {
            simplex.sort();
        }
        chain[dim] = simplices;
        for q in (0..dim).rev() {
            let mut q_simplex_set: BTreeSet<Simplex> = BTreeSet::new();
            for s in &chain[q + 1] {
                for idx in 0..(q + 2) {
                    let mut s = s.clone();
                    s.remove(idx);
                    s.sort();
                    q_simplex_set.insert(s);
                }
            }
            let q_chain: Vec<Simplex> = q_simplex_set.into_iter().collect();
            chain[q] = q_chain;
        }
        SimplicialComplex { chain: chain }
    }
    fn to_string(&self) -> String {
        let mut s = String::new();
        for (q, simplices) in self.chain.iter().enumerate() {
            s = s + &format!("C_{}:{:?}\n", q, simplices);
        }
        s
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
    fn boundary_mat(&self, q: usize) -> Option<Mat> {
        if q + 1 >= self.chain.len() {
            return None;
        }
        let mut mat = Mat::from_element(self.chain[q].len(), self.chain[q + 1].len(), 0);
        for (j, s) in self.chain[q + 1].iter().enumerate() {
            let result: Vec<(Simplex, i64)> = SimplicialComplex::delta(&s);
            for (s, sign) in result {
                let i = self.to_idx(&s).expect("something wrong");
                mat[(i, j)] = sign;
            }
        }
        Some(mat)
    }
    fn homology_group(&self, q: usize) -> FinitelyGeneratedModule {
        if q == 0 {
            let mat = self.boundary_mat(q).unwrap();
            let v = smith(&mat);
            let rank = mat.nrows() - v.len();
            FinitelyGeneratedModule {
                torsion: Vec::new(),
                rank: rank,
            }
        } else if q == self.chain.len() - 1 {
            let mat = self.boundary_mat(q - 1).unwrap();
            let v = smith(&mat);
            let rank = mat.ncols() - v.len();
            FinitelyGeneratedModule {
                torsion: Vec::new(),
                rank: rank,
            }
        } else {
            let mat1 = self.boundary_mat(q - 1).unwrap();
            let mat2 = self.boundary_mat(q).unwrap();
            let (v1, v2) = (smith(&mat1), smith(&mat2));
            let rank = mat1.ncols() - v1.len() - v2.len();
            let mut torsion: Vec<i64> = Vec::new();
            for a in v2 {
                if a > 1 {
                    torsion.push(a);
                }
            }
            FinitelyGeneratedModule {
                torsion: torsion,
                rank: rank,
            }
        }
    }
}

struct FinitelyGeneratedModule {
    torsion: Vec<i64>,
    rank: usize,
}

impl FinitelyGeneratedModule {
    fn to_string(&self) -> String {
        let mut torsion_str: Vec<String> =
            self.torsion.iter().map(|a| format!("Z/{}Z", a)).collect();
        if self.torsion.len() == 0 && self.rank == 0 {
            return "0".to_string();
        }
        let z_str = match self.rank {
            1 => "Z".to_string(),
            _ => format!("Z^{}", self.rank),
        };
        if self.rank > 0 {
            torsion_str.push(z_str);
        }
        torsion_str.join("⊕")
    }
}

fn read_file(filename: &str) -> Vec<Simplex> {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut simplices: Vec<Simplex> = Vec::new();
    for line in reader.lines() {
        let simplex: Simplex = match line {
            Ok(line) => line.trim()
                .split(' ')
                .map(|s| s.parse::<usize>().unwrap())
                .collect(),
            Err(error) => {
                panic!("{}", error);
            }
        };
        simplices.push(simplex);
    }
    simplices
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        println!("{} [simplices file]", args[0]);
        return;
    }
    let c2 = read_file(&args[1]);
    let sc = SimplicialComplex::from_simplices(&c2);
    println!("{}", sc.to_string());
    let mut euler_char: i64 = 0;
    for q in 0..sc.chain.len() {
        let h = sc.homology_group(q);
        println!("H_{} = {}", q, h.to_string());
        if q % 2 == 0 {
            euler_char = euler_char + (h.rank as i64);
        } else {
            euler_char = euler_char - (h.rank as i64);
        }
    }
    println!("Euler characteristic χ = {}", euler_char);
}
