# C' and C' - C on semi-perfect quantum codes

Every number below is computed on the 1909 semi-perfect codes only (q = 2, 3, 4, 5, 7, 8); normal codes are excluded throughout.


## 0. The exact identity

```
C' - C = [H_q2(delta_avg) - R] - [H_q2(delta_max) - R]
       =  H_q2(delta_avg) - H_q2(delta_max)
```
R cancels identically, so C' - C does not depend on X, on the rate, or
on n except through the two weight ratios. Verified numerically:
max |(C' - C) - (H_avg - H_max)| = 0.000e+00 over all codes.

Because delta_avg <= delta_max always, and H_q2 rises to a peak at
delta* = 1 - 1/q^2 then falls, three exact characterisations follow -
each confirmed on 100% of the relevant codes:

| condition | equivalent to | codes |
|---|---|---|
| C' - C = 0 | uniform weight spectrum (delta_avg = delta_max) | 1289 |
| C' - C > 0 | delta_max past the entropy peak delta* | 344 |
| C' = 0 | delta_max = 1 (full-support stabilizer rows) | 623 |

## 1. Structure by q

| q | codes | past peak | C'=0 | C'>0 | D=0 | D<0 | D>0 | min D | max D |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 384 | 89.3% | 33.3% | 65.4% | 38.5% | 13.5% | 47.9% | -0.09194 | +0.02797 |
| 3 | 249 | 67.9% | 39.8% | 60.2% | 50.6% | 28.9% | 20.5% | -0.01421 | +0.01410 |
| 4 | 236 | 76.3% | 41.9% | 56.8% | 55.1% | 24.2% | 20.8% | -0.01755 | +0.00692 |
| 5 | 245 | 64.5% | 40.4% | 59.2% | 70.2% | 10.6% | 19.2% | -0.00606 | +0.00405 |
| 7 | 353 | 41.4% | 28.0% | 66.6% | 86.4% | 9.9% | 3.7% | -0.00298 | +0.00067 |
| 8 | 442 | 30.3% | 22.4% | 74.0% | 92.3% | 7.7% | 0.0% | -0.00200 | +0.00000 |

Semi-perfect codes are dominated by uniform-weight stabilizer matrices, and increasingly so with q: from 39% at q=2 to 92% at q=8. For those codes C' and C coincide exactly.


## 2. Power-law fits  y = alpha n^beta

Fitted on strictly positive values only (the exact zeros above cannot
enter a log-log fit); `usable` is false where the fitted curve is worse
than a constant in linear space.

| series | q | target | alpha | beta | R2 log-log | R2 linear | points | usable |
|---|---|---|---|---|---|---|---|---|
| C' | 2 | max | 0.4351 | -0.8414 | 0.8020 | 0.8446 | 243 | yes |
| C' | 2 | mean | 0.3825 | -0.8156 | 0.7986 | 0.8515 | 243 | yes |
| C' | 3 | max | 0.1508 | -0.7018 | 0.5789 | 0.6551 | 93 | yes |
| C' | 3 | mean | 0.1420 | -0.7013 | 0.6115 | 0.6679 | 93 | yes |
| C' | 4 | max | 0.1415 | -0.9101 | 0.4047 | 0.5493 | 90 | yes |
| C' | 4 | mean | 0.2103 | -1.0606 | 0.5235 | 0.6771 | 90 | yes |
| C' | 5 | max | 0.7659 | -1.3677 | 0.7575 | 0.7466 | 95 | yes |
| C' | 5 | mean | 0.5943 | -1.3091 | 0.7582 | 0.8563 | 95 | yes |
| C' | 7 | max | 8.5870 | -2.1435 | 0.6962 | -16.2940 | 77 | NO |
| C' | 7 | mean | 6.4034 | -2.0896 | 0.7077 | -11.5732 | 77 | NO |
| C' | 8 | max | 6.0120 | -1.9782 | 0.5674 | -16.7255 | 81 | NO |
| C' | 8 | mean | 4.8289 | -1.9527 | 0.5862 | -12.5367 | 81 | NO |
| |C'-C| | 2 | max | 0.1107 | -1.0978 | 0.1505 | 0.0333 | 230 | yes |
| |C'-C| | 2 | mean | 0.0949 | -1.0675 | 0.1456 | 0.0632 | 230 | yes |
| |C'-C| | 3 | max | 0.2695 | -1.3876 | 0.2304 | -0.0308 | 80 | NO |
| |C'-C| | 3 | mean | 0.2840 | -1.4414 | 0.2617 | 0.0333 | 80 | yes |
| |C'-C| | 4 | max | 0.0029 | -0.1680 | 0.0066 | -0.0690 | 74 | NO |
| |C'-C| | 4 | mean | 0.0067 | -0.4198 | 0.0438 | -0.0117 | 74 | NO |
| |C'-C| | 5 | max | 0.0184 | -0.7975 | 0.0769 | -0.0313 | 73 | NO |
| |C'-C| | 5 | mean | 0.0184 | -0.7975 | 0.0769 | -0.0313 | 73 | NO |
| |C'-C| | 7 | max | 12178.2277 | -3.9151 | 0.4114 | 0.6066 | 48 | yes |
| |C'-C| | 7 | mean | 12178.2277 | -3.9151 | 0.4114 | 0.6066 | 48 | yes |
| |C'-C| | 8 | max | 11444.4376 | -3.7330 | 0.4648 | 0.3184 | 34 | yes |
| |C'-C| | 8 | mean | 11444.4376 | -3.7330 | 0.4648 | 0.3184 | 34 | yes |

## 3. Envelope

| series | q | (a ln n + b)/n | R2 | alpha n^beta | R2 lin | a*ln q2 |
|---|---|---|---|---|---|---|
| C' | 2 | (0.0547 ln n +0.6642)/n | 0.8981 | 0.4351 n^-0.8414 | 0.8446 | 0.0758 |
| C' | 3 | (-0.0020 ln n +0.4190)/n | 0.8387 | 0.1508 n^-0.7018 | 0.6551 | -0.0044 |
| C' | 4 | (-0.0470 ln n +0.4347)/n | 0.8727 | 0.1415 n^-0.9101 | 0.5493 | -0.1304 |
| C' | 5 | (0.0077 ln n +0.3168)/n | 0.8844 | 0.7659 n^-1.3677 | 0.7466 | 0.0248 |
| C' | 7 | (0.0669 ln n +0.1567)/n | 0.9127 | 8.5870 n^-2.1435 | -16.2940 | 0.2604 |
| C' | 8 | (0.0760 ln n +0.1207)/n | 0.9375 | 6.0120 n^-1.9782 | -16.7255 | 0.3159 |
| |C'-C| | 2 | (-0.0001 ln n +0.3518)/n | 0.2241 | 0.1107 n^-1.0978 | 0.0333 | -0.0002 |
| |C'-C| | 3 | (0.0513 ln n -0.0549)/n | 0.1448 | 0.2695 n^-1.3876 | -0.0308 | 0.1127 |
| |C'-C| | 4 | (-0.0019 ln n +0.1170)/n | 0.1433 | 0.0029 n^-0.1680 | -0.0690 | -0.0052 |
| |C'-C| | 5 | (0.0043 ln n +0.0462)/n | 0.1193 | 0.0184 n^-0.7975 | -0.0313 | 0.0139 |
| |C'-C| | 7 | (-0.1547 ln n +0.7269)/n | 0.6518 | 12178.2277 n^-3.9151 | 0.6066 | -0.6022 |
| |C'-C| | 8 | (-0.1476 ln n +0.7284)/n | 0.3834 | 11444.4376 n^-3.7330 | 0.3184 | -0.6139 |

## 4. Which variables explain C' and C' - C

| target | regime | q | features | R2 | RMSE | codes |
|---|---|---|---|---|---|---|
| C' | all semi-perfect | 2 | n only | 0.3538 | 0.01412 | 384 |
| C' | all semi-perfect | 2 | n, k | 0.7442 | 0.00889 | 384 |
| C' | all semi-perfect | 2 | n, k + delta | 0.9059 | 0.00539 | 384 |
| C' | all semi-perfect | 3 | n only | 0.2725 | 0.01218 | 249 |
| C' | all semi-perfect | 3 | n, k | 0.8712 | 0.00513 | 249 |
| C' | all semi-perfect | 3 | n, k + delta | 0.9156 | 0.00415 | 249 |
| C' | all semi-perfect | 4 | n only | 0.4006 | 0.01020 | 236 |
| C' | all semi-perfect | 4 | n, k | 0.8979 | 0.00421 | 236 |
| C' | all semi-perfect | 4 | n, k + delta | 0.9430 | 0.00315 | 236 |
| C' | all semi-perfect | 5 | n only | 0.5366 | 0.00945 | 245 |
| C' | all semi-perfect | 5 | n, k | 0.9476 | 0.00318 | 245 |
| C' | all semi-perfect | 5 | n, k + delta | 0.9517 | 0.00305 | 245 |
| C' | all semi-perfect | 7 | n only | 0.5338 | 0.00690 | 353 |
| C' | all semi-perfect | 7 | n, k | 0.9268 | 0.00273 | 353 |
| C' | all semi-perfect | 7 | n, k + delta | 0.9342 | 0.00259 | 353 |
| C' | all semi-perfect | 8 | n only | 0.4986 | 0.00588 | 442 |
| C' | all semi-perfect | 8 | n, k | 0.9244 | 0.00228 | 442 |
| C' | all semi-perfect | 8 | n, k + delta | 0.9318 | 0.00217 | 442 |
| C' | all semi-perfect | pooled | n only | 0.3374 | 0.01064 | 1909 |
| C' | all semi-perfect | pooled | n, k | 0.7102 | 0.00703 | 1909 |
| C' | all semi-perfect | pooled | n, k + delta | 0.7867 | 0.00603 | 1909 |
| C' | all semi-perfect | pooled | n, q | 0.3891 | 0.01021 | 1909 |
| C' | all semi-perfect | pooled | n, k, q | 0.7678 | 0.00630 | 1909 |
| C' | all semi-perfect | pooled | n, k, q + delta | 0.8230 | 0.00550 | 1909 |
| C' | below entropy peak | 2 | n only | 0.6238 | 0.02373 | 41 |
| C' | below entropy peak | 2 | n, k | 0.8775 | 0.01354 | 41 |
| C' | below entropy peak | 2 | n, k + delta | 0.9877 | 0.00430 | 41 |
| C' | below entropy peak | 3 | n only | 0.9084 | 0.00624 | 80 |
| C' | below entropy peak | 3 | n, k | 0.9498 | 0.00462 | 80 |
| C' | below entropy peak | 3 | n, k + delta | 0.9823 | 0.00274 | 80 |
| C' | below entropy peak | 4 | n only | 0.9204 | 0.00585 | 56 |
| C' | below entropy peak | 4 | n, k | 0.9461 | 0.00481 | 56 |
| C' | below entropy peak | 4 | n, k + delta | 0.9942 | 0.00158 | 56 |
| C' | below entropy peak | 5 | n only | 0.9151 | 0.00403 | 87 |
| C' | below entropy peak | 5 | n, k | 0.9887 | 0.00147 | 87 |
| C' | below entropy peak | 5 | n, k + delta | 0.9961 | 0.00087 | 87 |
| C' | below entropy peak | 7 | n only | 0.9201 | 0.00257 | 207 |
| C' | below entropy peak | 7 | n, k | 0.9942 | 0.00069 | 207 |
| C' | below entropy peak | 7 | n, k + delta | 0.9953 | 0.00062 | 207 |
| C' | below entropy peak | 8 | n only | 0.9298 | 0.00204 | 308 |
| C' | below entropy peak | 8 | n, k | 0.9962 | 0.00047 | 308 |
| C' | below entropy peak | 8 | n, k + delta | 0.9962 | 0.00047 | 308 |
| C' | below entropy peak | pooled | n only | 0.6638 | 0.00913 | 779 |
| C' | below entropy peak | pooled | n, k | 0.7165 | 0.00839 | 779 |
| C' | below entropy peak | pooled | n, k + delta | 0.8081 | 0.00690 | 779 |
| C' | below entropy peak | pooled | n, q | 0.8157 | 0.00676 | 779 |
| C' | below entropy peak | pooled | n, k, q | 0.8797 | 0.00546 | 779 |
| C' | below entropy peak | pooled | n, k, q + delta | 0.9708 | 0.00269 | 779 |
| C' | past entropy peak | 2 | n only | 0.2387 | 0.00900 | 343 |
| C' | past entropy peak | 2 | n, k | 0.8477 | 0.00403 | 343 |
| C' | past entropy peak | 2 | n, k + delta | 0.9124 | 0.00305 | 343 |
| C' | past entropy peak | 3 | n only | 0.0153 | 0.00579 | 169 |
| C' | past entropy peak | 3 | n, k | 0.7394 | 0.00298 | 169 |
| C' | past entropy peak | 3 | n, k + delta | 0.7569 | 0.00288 | 169 |
| C' | past entropy peak | 4 | n only | 0.0073 | 0.00212 | 180 |
| C' | past entropy peak | 4 | n, k | 0.3121 | 0.00176 | 180 |
| C' | past entropy peak | 4 | n, k + delta | 0.6072 | 0.00133 | 180 |
| C' | past entropy peak | 5 | n only | 0.0360 | 0.00148 | 158 |
| C' | past entropy peak | 5 | n, k | 0.5633 | 0.00099 | 158 |
| C' | past entropy peak | 5 | n, k + delta | 0.6437 | 0.00090 | 158 |
| C' | past entropy peak | 7 | n only | 0.0025 | 0.00052 | 146 |
| C' | past entropy peak | 7 | n, k | 0.0265 | 0.00051 | 146 |
| C' | past entropy peak | 7 | n, k + delta | 0.7757 | 0.00025 | 146 |
| C' | past entropy peak | 8 | n only | 0.0001 | 0.00038 | 134 |
| C' | past entropy peak | 8 | n, k | 0.0051 | 0.00038 | 134 |
| C' | past entropy peak | 8 | n, k + delta | 0.9396 | 0.00009 | 134 |
| C' | past entropy peak | pooled | n only | 0.0075 | 0.00669 | 1130 |
| C' | past entropy peak | pooled | n, k | 0.6690 | 0.00386 | 1130 |
| C' | past entropy peak | pooled | n, k + delta | 0.8027 | 0.00298 | 1130 |
| C' | past entropy peak | pooled | n, q | 0.2980 | 0.00563 | 1130 |
| C' | past entropy peak | pooled | n, k, q | 0.7183 | 0.00356 | 1130 |
| C' | past entropy peak | pooled | n, k, q + delta | 0.8135 | 0.00290 | 1130 |
| C'-C | all semi-perfect | 2 | n only | 0.0405 | 0.00977 | 384 |
| C'-C | all semi-perfect | 2 | n, k | 0.6248 | 0.00611 | 384 |
| C'-C | all semi-perfect | 2 | n, k + delta | 0.7285 | 0.00520 | 384 |
| C'-C | all semi-perfect | 3 | n only | 0.0187 | 0.00287 | 249 |
| C'-C | all semi-perfect | 3 | n, k | 0.0518 | 0.00282 | 249 |
| C'-C | all semi-perfect | 3 | n, k + delta | 0.1161 | 0.00272 | 249 |
| C'-C | all semi-perfect | 4 | n only | 0.0065 | 0.00221 | 236 |
| C'-C | all semi-perfect | 4 | n, k | 0.0450 | 0.00217 | 236 |
| C'-C | all semi-perfect | 4 | n, k + delta | 0.0556 | 0.00216 | 236 |
| C'-C | all semi-perfect | 5 | n only | 0.0394 | 0.00084 | 245 |
| C'-C | all semi-perfect | 5 | n, k | 0.0482 | 0.00083 | 245 |
| C'-C | all semi-perfect | 5 | n, k + delta | 0.1440 | 0.00079 | 245 |
| C'-C | all semi-perfect | 7 | n only | 0.0393 | 0.00041 | 353 |
| C'-C | all semi-perfect | 7 | n, k | 0.0424 | 0.00041 | 353 |
| C'-C | all semi-perfect | 7 | n, k + delta | 0.4736 | 0.00030 | 353 |
| C'-C | all semi-perfect | 8 | n only | 0.0630 | 0.00028 | 442 |
| C'-C | all semi-perfect | 8 | n, k | 0.0833 | 0.00028 | 442 |
| C'-C | all semi-perfect | 8 | n, k + delta | 0.8552 | 0.00011 | 442 |
| C'-C | all semi-perfect | pooled | n only | 0.0065 | 0.00466 | 1909 |
| C'-C | all semi-perfect | pooled | n, k | 0.1443 | 0.00433 | 1909 |
| C'-C | all semi-perfect | pooled | n, k + delta | 0.1891 | 0.00421 | 1909 |
| C'-C | all semi-perfect | pooled | n, q | 0.0262 | 0.00462 | 1909 |
| C'-C | all semi-perfect | pooled | n, k, q | 0.2170 | 0.00414 | 1909 |
| C'-C | all semi-perfect | pooled | n, k, q + delta | 0.2130 | 0.00415 | 1909 |
| C'-C | below entropy peak | 2 | n only | 0.1931 | 0.02338 | 41 |
| C'-C | below entropy peak | 2 | n, k | 0.8714 | 0.00934 | 41 |
| C'-C | below entropy peak | 2 | n, k + delta | 0.9804 | 0.00364 | 41 |
| C'-C | below entropy peak | 3 | n only | 0.2245 | 0.00229 | 80 |
| C'-C | below entropy peak | 3 | n, k | 0.3071 | 0.00217 | 80 |
| C'-C | below entropy peak | 3 | n, k + delta | 0.8689 | 0.00094 | 80 |
| C'-C | below entropy peak | 4 | n only | 0.0546 | 0.00307 | 56 |
| C'-C | below entropy peak | 4 | n, k | 0.0743 | 0.00304 | 56 |
| C'-C | below entropy peak | 4 | n, k + delta | 0.8643 | 0.00116 | 56 |
| C'-C | below entropy peak | 5 | n only | 0.0780 | 0.00071 | 87 |
| C'-C | below entropy peak | 5 | n, k | 0.0882 | 0.00070 | 87 |
| C'-C | below entropy peak | 5 | n, k + delta | 0.8844 | 0.00025 | 87 |
| C'-C | below entropy peak | 7 | n only | 0.0217 | 0.00015 | 207 |
| C'-C | below entropy peak | 7 | n, k | 0.0323 | 0.00015 | 207 |
| C'-C | below entropy peak | 7 | n, k + delta | 1.0000 | 0.00000 | 207 |
| C'-C | below entropy peak | 8 | n only | nan | 0.00000 | 308 |
| C'-C | below entropy peak | 8 | n, k | nan | 0.00000 | 308 |
| C'-C | below entropy peak | 8 | n, k + delta | nan | 0.00000 | 308 |
| C'-C | below entropy peak | pooled | n only | 0.0056 | 0.00670 | 779 |
| C'-C | below entropy peak | pooled | n, k | 0.2796 | 0.00571 | 779 |
| C'-C | below entropy peak | pooled | n, k + delta | 0.7856 | 0.00311 | 779 |
| C'-C | below entropy peak | pooled | n, q | 0.2478 | 0.00583 | 779 |
| C'-C | below entropy peak | pooled | n, k, q | 0.5774 | 0.00437 | 779 |
| C'-C | below entropy peak | pooled | n, k, q + delta | 0.8017 | 0.00299 | 779 |
| C'-C | past entropy peak | 2 | n only | 0.0419 | 0.00313 | 343 |
| C'-C | past entropy peak | 2 | n, k | 0.1613 | 0.00293 | 343 |
| C'-C | past entropy peak | 2 | n, k + delta | 0.4074 | 0.00247 | 343 |
| C'-C | past entropy peak | 3 | n only | 0.0078 | 0.00268 | 169 |
| C'-C | past entropy peak | 3 | n, k | 0.0591 | 0.00261 | 169 |
| C'-C | past entropy peak | 3 | n, k + delta | 0.3520 | 0.00217 | 169 |
| C'-C | past entropy peak | 4 | n only | 0.0171 | 0.00163 | 180 |
| C'-C | past entropy peak | 4 | n, k | 0.0373 | 0.00162 | 180 |
| C'-C | past entropy peak | 4 | n, k + delta | 0.1308 | 0.00154 | 180 |
| C'-C | past entropy peak | 5 | n only | 0.0461 | 0.00085 | 158 |
| C'-C | past entropy peak | 5 | n, k | 0.1494 | 0.00080 | 158 |
| C'-C | past entropy peak | 5 | n, k + delta | 0.2597 | 0.00075 | 158 |
| C'-C | past entropy peak | 7 | n only | 0.0187 | 0.00060 | 146 |
| C'-C | past entropy peak | 7 | n, k | 0.5216 | 0.00042 | 146 |
| C'-C | past entropy peak | 7 | n, k + delta | 0.5400 | 0.00041 | 146 |
| C'-C | past entropy peak | 8 | n only | 0.0894 | 0.00047 | 134 |
| C'-C | past entropy peak | 8 | n, k | 0.7517 | 0.00025 | 134 |
| C'-C | past entropy peak | 8 | n, k + delta | 0.8394 | 0.00020 | 134 |
| C'-C | past entropy peak | pooled | n only | 0.0010 | 0.00224 | 1130 |
| C'-C | past entropy peak | pooled | n, k | 0.0890 | 0.00214 | 1130 |
| C'-C | past entropy peak | pooled | n, k + delta | 0.2281 | 0.00197 | 1130 |
| C'-C | past entropy peak | pooled | n, q | 0.0630 | 0.00217 | 1130 |
| C'-C | past entropy peak | pooled | n, k, q | 0.1061 | 0.00212 | 1130 |
| C'-C | past entropy peak | pooled | n, k, q + delta | 0.2560 | 0.00194 | 1130 |
