Sagemath reference implementation of the Arora-Ge algorithm
 for breaking NTRU with multiple keys.

To load the package in the REPL for interactive use do 
`attach("arora_ge_ntru.sage")`.

# Usage

To generate and solve a multiple-key NTRU instance in the degree `n` ring 
modulo `q` with `m` samples do:
```
$ sage arora_ge_ntru.sage n q m [-c COEFFS] [-r RING]
```
where `r = 1` for `Z[x]/(x^n - 1)` or `r = 2` for `Z[x]/(x^n + 1)`. Set `c = 2`
for binary support or `c = 3` for ternary. 

For example, this is the same parameters as the first entry of table 1:
```
$ sage arora_ge_ntru.sage 32 769 18 -c 2 -r 1

Building linear system... 6.367090419167653 seconds.
System dimensions: 576 x 560.

Computing right kernel... 0.05066058994270861 seconds.
Kernel dimensions: 32 x 560

Doing kernel reduction...
Computing kernel of 32 x 32 submatrix...
Initial kernel rank: 9.

Computing kernel of 32 x 62 submatrix...
Reduced kernel rank: 4.

Computing kernel of 32 x 91 submatrix...
Reduced kernel rank: 1.

SUCCESS.
Kernel reduction time: 0.008852346800267696 seconds.

Recovered gg with gg = g*x^19.
```


