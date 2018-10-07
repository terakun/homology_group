# homology_group
calculate simplicial complex homology group
### How to use
```
cargo run [simplices file]
```
### Example
projective plane:

![projective_plane](https://i.imgur.com/pxpFO22.png)

Corresponding file:
`projectiveplane.dat`
```
1 2 3 
2 0 3 
1 4 3 
1 0 4 
0 2 4 
4 5 3 
3 5 0 
0 5 1 
4 2 5 
5 2 1 
```

result:
```
$ cargo run projectiveplane.dat 
C_0:[[0], [1], [2], [3], [4], [5]]
C_1:[[0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [1, 2], [1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]]
C_2:[[1, 2, 3], [0, 2, 3], [1, 3, 4], [0, 1, 4], [0, 2, 4], [3, 4, 5], [0, 3, 5], [0, 1, 5], [2, 4, 5], [1, 2, 5]]

H_0 = Z
H_1 = Z/2Z
H_2 = 0
Euler characteristic χ = 1
```
