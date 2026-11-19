<!-- https://docs.dseams.info/explanation/clavis-trace.md -->

# ClavisCage TRACE {#claviscage-trace}

ClavisCage 0.3.18 runs TRACE (Hsu and Lin, J. Chem. Theory Comput. 2025,
[10.1021/acs.jctc.5c01459](https://doi.org/10.1021/acs.jctc.5c01459))
on the same hydrate cells as d-SEAMS 2.14.0. d-SEAMS counts closed cages
of census $5^{12}$, $5^{12}6^{2}$ and $5^{12}6^{4}$ on the oxygen
neighbour graph. TRACE counts cages on the hydrogen-bond graph. The two
constructions are [TRACE and signature cages](trace-assembly.md).
The commands are the how-to
[Cages and TRACE](../user/howto/cages-trace.md).

`analyze_waters` takes no thread count, so an eight-thread row times the
same serial TRACE call as the one-thread row. d-SEAMS uses the width the
row names. A row also names the version pair, the molecule count, the
repeat count `n` and the core count. The agreement text is the cage
census: cages found by both codes, cages found only by d-SEAMS, and
cages found only by TRACE.

## ClavisCage TRACE {#id1}

The tables compare d-SEAMS with the named analysis packages on the same inputs.

One row per cell. The agreement column is the cage census. The signature-by-signature rows are in the block below.

On the 46-molecule sI cell the census drops 7 cages that only d-SEAMS counts, and TRACE counts none that d-SEAMS does not.
512-0 has 18 of 30 bonds across the periodic cell; the six 51262 cages have 10 to 15. The 5^12 TRACE returns has none. The 2x2x2 tile matches, so d-SEAMS is right.
On 110160 molecules, one thread, d-SEAMS takes 955.8314 ms and TRACE takes 61962.4841 ms. At eight threads d-SEAMS takes 303.0753 ms and TRACE takes 61956.9697 ms. TRACE is serial: analyze_waters has no thread argument and Cargo.toml leaves rayon commented out.

| version | molecules | threads | cores | n | d-SEAMS ms | TRACE ms | ratio | agreement |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.14.0/0.3.18 | 46 | 1 | 8 | 8 | 0.5796 | 50.3899 | 86.9321 | 512 both 1 d-SEAMS-only 1 TRACE-only 0; 51262 both 0 d-SEAMS-only 6 TRACE-only 0; 51264 both 0 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 46 | 8 | 8 | 8 | 0.5947 | 50.2893 | 84.5632 | 512 both 1 d-SEAMS-only 1 TRACE-only 0; 51262 both 0 d-SEAMS-only 6 TRACE-only 0; 51264 both 0 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 368 | 1 | 8 | 8 | 3.7222 | 297.6070 | 79.9538 | 512 both 16 d-SEAMS-only 0 TRACE-only 0; 51262 both 48 d-SEAMS-only 0 TRACE-only 0; 51264 both 0 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 368 | 8 | 8 | 8 | 2.9243 | 297.9639 | 101.8907 | 512 both 16 d-SEAMS-only 0 TRACE-only 0; 51262 both 48 d-SEAMS-only 0 TRACE-only 0; 51264 both 0 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 136 | 1 | 8 | 8 | 1.4924 | 68.3007 | 45.7666 | 512 both 16 d-SEAMS-only 0 TRACE-only 0; 51262 both 0 d-SEAMS-only 0 TRACE-only 0; 51264 both 8 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 136 | 8 | 8 | 8 | 1.5655 | 68.2571 | 43.6016 | 512 both 16 d-SEAMS-only 0 TRACE-only 0; 51262 both 0 d-SEAMS-only 0 TRACE-only 0; 51264 both 8 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 1088 | 1 | 8 | 8 | 10.7150 | 697.6882 | 65.1132 | 512 both 128 d-SEAMS-only 0 TRACE-only 0; 51262 both 0 d-SEAMS-only 0 TRACE-only 0; 51264 both 64 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 1088 | 8 | 8 | 8 | 4.7284 | 697.7931 | 147.5752 | 512 both 128 d-SEAMS-only 0 TRACE-only 0; 51262 both 0 d-SEAMS-only 0 TRACE-only 0; 51264 both 64 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 110160 | 1 | 8 | 8 | 955.8314 | 61962.4841 | 64.8257 | 512 both 12960 d-SEAMS-only 0 TRACE-only 0; 51262 both 0 d-SEAMS-only 0 TRACE-only 0; 51264 both 6480 d-SEAMS-only 0 TRACE-only 0 |
| 2.14.0/0.3.18 | 110160 | 8 | 8 | 8 | 303.0753 | 61956.9697 | 204.4277 | 512 both 12960 d-SEAMS-only 0 TRACE-only 0; 51262 both 0 d-SEAMS-only 0 TRACE-only 0; 51264 both 6480 d-SEAMS-only 0 TRACE-only 0 |

| version | molecules | threads | cores | n | signature | both | d-SEAMS only | TRACE only |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2.14.0/0.3.18 | 46 | 1 | 8 | 8 | 512 | 1 | 1 | 0 |
| 2.14.0/0.3.18 | 46 | 1 | 8 | 8 | 51262 | 0 | 6 | 0 |
| 2.14.0/0.3.18 | 46 | 1 | 8 | 8 | 51264 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 46 | 8 | 8 | 8 | 512 | 1 | 1 | 0 |
| 2.14.0/0.3.18 | 46 | 8 | 8 | 8 | 51262 | 0 | 6 | 0 |
| 2.14.0/0.3.18 | 46 | 8 | 8 | 8 | 51264 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 368 | 1 | 8 | 8 | 512 | 16 | 0 | 0 |
| 2.14.0/0.3.18 | 368 | 1 | 8 | 8 | 51262 | 48 | 0 | 0 |
| 2.14.0/0.3.18 | 368 | 1 | 8 | 8 | 51264 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 368 | 8 | 8 | 8 | 512 | 16 | 0 | 0 |
| 2.14.0/0.3.18 | 368 | 8 | 8 | 8 | 51262 | 48 | 0 | 0 |
| 2.14.0/0.3.18 | 368 | 8 | 8 | 8 | 51264 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 136 | 1 | 8 | 8 | 512 | 16 | 0 | 0 |
| 2.14.0/0.3.18 | 136 | 1 | 8 | 8 | 51262 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 136 | 1 | 8 | 8 | 51264 | 8 | 0 | 0 |
| 2.14.0/0.3.18 | 136 | 8 | 8 | 8 | 512 | 16 | 0 | 0 |
| 2.14.0/0.3.18 | 136 | 8 | 8 | 8 | 51262 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 136 | 8 | 8 | 8 | 51264 | 8 | 0 | 0 |
| 2.14.0/0.3.18 | 1088 | 1 | 8 | 8 | 512 | 128 | 0 | 0 |
| 2.14.0/0.3.18 | 1088 | 1 | 8 | 8 | 51262 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 1088 | 1 | 8 | 8 | 51264 | 64 | 0 | 0 |
| 2.14.0/0.3.18 | 1088 | 8 | 8 | 8 | 512 | 128 | 0 | 0 |
| 2.14.0/0.3.18 | 1088 | 8 | 8 | 8 | 51262 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 1088 | 8 | 8 | 8 | 51264 | 64 | 0 | 0 |
| 2.14.0/0.3.18 | 110160 | 1 | 8 | 8 | 512 | 12960 | 0 | 0 |
| 2.14.0/0.3.18 | 110160 | 1 | 8 | 8 | 51262 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 110160 | 1 | 8 | 8 | 51264 | 6480 | 0 | 0 |
| 2.14.0/0.3.18 | 110160 | 8 | 8 | 8 | 512 | 12960 | 0 | 0 |
| 2.14.0/0.3.18 | 110160 | 8 | 8 | 8 | 51262 | 0 | 0 | 0 |
| 2.14.0/0.3.18 | 110160 | 8 | 8 | 8 | 51264 | 6480 | 0 | 0 |
