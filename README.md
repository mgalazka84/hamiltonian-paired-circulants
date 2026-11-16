# Hamiltonian-paired circulant graphs of valency four

Computational supplement to the manuscript by Marek Gałązka and Hanna Wdowicka.

A two-jump circulant has vertices labelled from `0` to `n-1`. From a vertex `x`, one can move by `a`, `-a`, `b`, or `-b`, with labels taken modulo `n`. A Hamilton cycle visits every vertex once before returning to its start. In the graphs studied here, its unused edges form one or more cycles. The graph is Hamiltonian-paired when it has a Hamilton cycle and every such choice leaves a single cycle.

`verify.py` reproduces the finite verification in Section 6 and the example in Figure 3. It uses Python 3.10 or later and only the standard library. Run the commands below from this directory.

## Reproduce the counts in Section 6

```sh
python verify.py census
```

```text
Through order 80: 16824 parameter triples, 87 in the listed families, 16737 negative certificates checked.
Direct search through order 26: 508 triples, 9 paired, 499 non-paired.
```

The construction checks every triple with `5 ≤ n ≤ 80`, `1 ≤ a < b < n/2` and `gcd(n,a,b) = 1`. For each triple outside the classified families, it constructs a Hamilton cycle and checks that the unused edges form at least two cycles. The second line comes from an independent search through Hamilton cycles for `n ≤ 26`. The counts distinguish jump presentations, so different presentations can describe isomorphic graphs.

## Reproduce Figure 3

```sh
python verify.py witness 13 1 5
```

```text
Square steps: (1, 5) Origins: (0, 4)
Hamilton cycle: (0, 1, 9, 10, 2, 7, 12, 4, 5, 6, 11, 3, 8)
Complement cycles: ((0, 5, 10, 11, 12), (1, 2, 3, 4, 9, 8, 7, 6))
```

Each vertex list closes by returning to its first entry. The unused edges in this example form a cycle of length five and a cycle of length eight. To reproduce the example with unequal initial cycle counts in Section 3, run `python verify.py witness 12 3 4`.

For a shorter check, run `python verify.py census --max-order 20 --exhaustive-order 12`.
