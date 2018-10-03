extern crate nalgebra as na;
use na::DMatrix;
type Mat = DMatrix<i64>;
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
    let mat1 = sc.boundary_mat(0);
    let mat2 = sc.boundary_mat(1);
    println!("d1 = {}", mat1);
    println!("d2 = {}", mat2);
    let mat3 = mat1 * mat2;
    println!("d1 * d2 = {}", mat3);
}
