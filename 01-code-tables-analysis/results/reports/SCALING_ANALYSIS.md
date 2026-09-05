# Scaling of the correction factors C and C'

Fitted from 34775 codes across q = 2, 3, 4, 5, 7, 8 (see `analysis/` for the full tables and plots).


## Headline

**C is not a power law in n.** Expanding the definitions with Stirling's
formula gives, for a code whose largest stabilizer weight is
delta = delta_max,

```
    X       = sum_w  binom(n,w) (q^2-1)^w      >=  binom(n,wmax)(q^2-1)^wmax
    R       = log_{q^2}(X) / n
    log_{q^2}[ binom(n,w)(q^2-1)^w ] = n H_{q^2}(w/n)
                                       - (1/2) log_{q^2}(2 pi n d(1-d)) + O(1/n)
=>  C = H_{q^2}(delta) - R
      = log_{q^2}( 2 pi n delta(1-delta) ) / (2n)      <-- leading term
        - (1/n) log_{q^2}( X / X_dominant )            <-- correction
```
So the natural law is **C = Theta(log n / n)**, i.e. C = o(n^(e-1)) for
every e > 0 - not alpha*n^beta. A power law fitted over 2 <= n <= 256
nevertheless looks excellent (R2 ~ 0.99 in log-log) because log n is
nearly constant over two decades; the fitted beta ~ -0.75 is the log
factor masquerading as a shifted exponent, not a real exponent.

The behaviour splits at the entropy peak delta* = 1 - 1/q^2:

* **delta_max <= delta***: the largest weight is also the dominant term,
  C > 0 and C = Theta(log n / n).
* **delta_max > delta***: some smaller weight dominates X, so
  H_{q^2}(delta_max) < H_{q^2}(delta_dominant) and C < 0, of order 1 and
  essentially independent of n. Every single negative-C code in the whole
  data set (100%, all q) lies in this regime.


## 1. Power-law fits  y = alpha * n^beta

Fitted by least squares in log-log space on the per-length upper envelope `max_k y(n,k)` and on the per-length mean. Reported because the question asks for them; see the Headline for why the functional form is really log(n)/n. `R2 (log-log)` is the quality of the fit that was actually minimised; `R2 (linear)` is how well the same curve explains the untransformed values, and it is where the poor cases show themselves (a negative value means the curve is worse than a constant).


**C (all codes)**

| q | target | alpha | beta | R2 (log-log) | R2 (linear) | points |
|---|---|---|---|---|---|---|
| 2 | max_k C(n,k) | 0.6138 | -0.7752 | 0.9927 | 0.9092 | 255 |
| 2 | mean_k C(n,k) | 0.5671 | -0.8779 | 0.9094 | 0.0411 | 255 |
| 3 | max_k C(n,k) | 0.3141 | -0.7145 | 0.9939 | 0.9786 | 98 |
| 3 | mean_k C(n,k) | 0.1096 | -0.5426 | 0.8751 | 0.6394 | 98 |
| 4 | max_k C(n,k) | 0.2707 | -0.7443 | 0.9879 | 0.9555 | 97 |
| 4 | mean_k C(n,k) | 0.1113 | -0.6139 | 0.8519 | 0.6127 | 97 |
| 5 | max_k C(n,k) | 0.2652 | -0.7836 | 0.9807 | 0.9301 | 97 |
| 5 | mean_k C(n,k) | 0.1242 | -0.6886 | 0.8245 | 0.4090 | 97 |
| 7 | max_k C(n,k) | 0.1903 | -0.7425 | 0.9869 | 0.9364 | 99 |
| 7 | mean_k C(n,k) | 0.1365 | -0.7290 | 0.9229 | 0.5698 | 99 |
| 8 | max_k C(n,k) | 0.1815 | -0.7488 | 0.9842 | 0.9283 | 99 |
| 8 | mean_k C(n,k) | 0.1014 | -0.6801 | 0.9060 | 0.6306 | 99 |

**C' (semi-perfect codes only)**

| q | target | alpha | beta | R2 (log-log) | R2 (linear) | points |
|---|---|---|---|---|---|---|
| 2 | max_k C'(n,k) | 0.4351 | -0.8414 | 0.8020 | 0.8446 | 243 |
| 2 | mean_k C'(n,k) | 0.3094 | -0.8426 | 0.7057 | 0.7776 | 240 |
| 3 | max_k C'(n,k) | 0.1508 | -0.7018 | 0.5789 | 0.6551 | 93 |
| 3 | mean_k C'(n,k) | 0.0764 | -0.6743 | 0.5102 | 0.5972 | 93 |
| 4 | max_k C'(n,k) | 0.1415 | -0.9101 | 0.4047 | 0.5493 | 90 |
| 4 | mean_k C'(n,k) | 0.1100 | -1.0298 | 0.4835 | 0.5891 | 89 |
| 5 | max_k C'(n,k) | 0.7659 | -1.3677 | 0.7575 | 0.7466 | 95 |
| 5 | mean_k C'(n,k) | 0.6013 | -1.4726 | 0.7563 | 0.5056 | 95 |
| 7 | max_k C'(n,k) | 8.5870 | -2.1435 | 0.6962 | -16.2940 | 77 |
| 7 | mean_k C'(n,k) | 7.5836 | -2.2503 | 0.6782 | -26.8643 | 77 |
| 8 | max_k C'(n,k) | 6.0120 | -1.9782 | 0.5674 | -16.7255 | 81 |
| 8 | mean_k C'(n,k) | 4.7398 | -2.0365 | 0.5370 | -25.9577 | 81 |

## 2. Envelope

Two envelopes are fitted to the per-length maximum `max_k y(n,k)`: the
log(n)/n form the derivation implies, and the requested power law.

| q | series | envelope (a ln n + b)/n | R2 | power law alpha n^beta | R2 |
|---|---|---|---|---|---|
| 2 | C | (0.4292 ln n +0.0395)/n | 0.9507 | 0.6138 n^-0.7752 | 0.9092 |
| 3 | C | (0.2612 ln n -0.0191)/n | 0.7995 | 0.3141 n^-0.7145 | 0.9786 |
| 4 | C | (0.1875 ln n -0.0268)/n | 0.5980 | 0.2707 n^-0.7443 | 0.9555 |
| 5 | C | (0.1632 ln n -0.0231)/n | 0.6115 | 0.2652 n^-0.7836 | 0.9301 |
| 7 | C | (0.1522 ln n +0.0131)/n | 0.9295 | 0.1903 n^-0.7425 | 0.9364 |
| 8 | C | (0.1424 ln n +0.0123)/n | 0.9293 | 0.1815 n^-0.7488 | 0.9283 |
| 2 | C' | (0.3167 ln n -0.1719)/n | 0.6079 | 0.4351 n^-0.8414 | 0.8446 |
| 3 | C' | (0.1803 ln n -0.0784)/n | 0.4439 | 0.1508 n^-0.7018 | 0.6551 |
| 4 | C' | (0.1342 ln n -0.0440)/n | 0.4856 | 0.1415 n^-0.9101 | 0.5493 |
| 5 | C' | (0.1547 ln n -0.0723)/n | 0.5908 | 0.7659 n^-1.3677 | 0.7466 |
| 7 | C' | (0.1514 ln n -0.0839)/n | 0.6666 | 8.5870 n^-2.1435 | -16.2940 |
| 8 | C' | (0.1515 ln n -0.0900)/n | 0.7137 | 6.0120 n^-1.9782 | -16.7255 |

### 2b. q-collapse of the fitted constants

The derivation puts a factor 1/log(q^2) in front of everything, so if
the law is real then alpha*ln(q^2) and a*ln(q^2) should be constant in
q. They are, to within a few percent - which is what licenses writing
C as a function of (n, q) rather than one fit per q.

| series | q | alpha | beta | alpha*ln(q2) | env a | a*ln(q2) | 1/(2 ln q2) |
|---|---|---|---|---|---|---|---|
| C | 2 | 0.6138 | -0.7752 | 0.8510 | 0.4292 | 0.5950 | 0.3607 |
| C | 3 | 0.3141 | -0.7145 | 0.6901 | 0.2612 | 0.5739 | 0.2276 |
| C | 4 | 0.2707 | -0.7443 | 0.7506 | 0.1875 | 0.5199 | 0.1803 |
| C | 5 | 0.2652 | -0.7836 | 0.8536 | 0.1632 | 0.5254 | 0.1553 |
| C | 7 | 0.1903 | -0.7425 | 0.7407 | 0.1522 | 0.5924 | 0.1285 |
| C | 8 | 0.1815 | -0.7488 | 0.7549 | 0.1424 | 0.5922 | 0.1202 |
| C' | 2 | 0.4351 | -0.8414 | 0.6031 | 0.3167 | 0.4391 | 0.3607 |
| C' | 3 | 0.1508 | -0.7018 | 0.3314 | 0.1803 | 0.3961 | 0.2276 |
| C' | 4 | 0.1415 | -0.9101 | 0.3923 | 0.1342 | 0.3721 | 0.1803 |
| C' | 5 | 0.7659 | -1.3677 | 2.4653 | 0.1547 | 0.4980 | 0.1553 |
| C' | 7 | 8.5870 | -2.1435 | 33.4190 | 0.1514 | 0.5893 | 0.1285 |
| C' | 8 | 6.0120 | -1.9782 | 25.0033 | 0.1515 | 0.6302 | 0.1202 |

### 2c. Parameter-free bound  C <= log_{q2}(pi n / 2) / (2n)

| q | codes | coverage (all) | coverage (below peak) | violations | max n |
|---|---|---|---|---|---|
| 2 | 17882 | 0.9968 | 0.9957 | 58 | 118 |
| 3 | 4511 | 0.9933 | 0.9924 | 30 | 65 |
| 4 | 3707 | 0.9962 | 0.9956 | 14 | 42 |
| 5 | 2891 | 0.9896 | 0.9881 | 30 | 44 |
| 7 | 2971 | 0.9906 | 0.9897 | 28 | 40 |
| 8 | 2813 | 0.9904 | 0.9894 | 27 | 42 |

## 3. Which variables explain C

Linear models (with intercept) on the features named; `delta` means the
single derived feature log_{q^2}(2 pi n d(1-d))/(2n), which uses the
code's actual weight spectrum and is therefore NOT a function of
(n,k,q). Within a fixed q, adding q as a predictor cannot help, so the
q-bearing feature sets are only scored on the pooled data.

Reading of the table:

* **Below the peak**, (n) alone already gives R2 ~ 0.83-0.91 per q, and
  the pooled (n,q) model reaches 0.89 - so C really can be plotted as a
  surface over (n,q). Adding k lifts it to ~0.95.
* **Past the peak**, (n) explains essentially nothing (R2 ~ 0.00-0.05)
  and even (n,k,q) only reaches 0.12 pooled. Which weight happens to
  dominate X is a property of the individual stabilizer matrix, not of
  the parameters, so no function of (n,k,q) can capture it. Supplying
  delta rescues it (0.37-0.76).

| regime | q | features | R2 | RMSE | codes |
|---|---|---|---|---|---|
| all codes | 2 | n only | 0.3806 | 0.01012 | 17882 |
| all codes | 2 | n, k | 0.5448 | 0.00867 | 17882 |
| all codes | 2 | n, delta (theory) | 0.0498 | 0.01253 | 17882 |
| all codes | 2 | n, k + delta | 0.6842 | 0.00722 | 17882 |
| all codes | 3 | n only | 0.2705 | 0.00844 | 4511 |
| all codes | 3 | n, k | 0.5087 | 0.00692 | 4511 |
| all codes | 3 | n, delta (theory) | 0.0601 | 0.00958 | 4511 |
| all codes | 3 | n, k + delta | 0.7825 | 0.00461 | 4511 |
| all codes | 4 | n only | 0.2764 | 0.00715 | 3707 |
| all codes | 4 | n, k | 0.6289 | 0.00512 | 3707 |
| all codes | 4 | n, delta (theory) | 0.0657 | 0.00813 | 3707 |
| all codes | 4 | n, k + delta | 0.8540 | 0.00321 | 3707 |
| all codes | 5 | n only | 0.3482 | 0.00619 | 2891 |
| all codes | 5 | n, k | 0.7372 | 0.00393 | 2891 |
| all codes | 5 | n, delta (theory) | 0.0634 | 0.00742 | 2891 |
| all codes | 5 | n, k + delta | 0.9221 | 0.00214 | 2891 |
| all codes | 7 | n only | 0.6540 | 0.00398 | 2971 |
| all codes | 7 | n, k | 0.8045 | 0.00299 | 2971 |
| all codes | 7 | n, delta (theory) | 0.0366 | 0.00664 | 2971 |
| all codes | 7 | n, k + delta | 0.9577 | 0.00139 | 2971 |
| all codes | 8 | n only | 0.4855 | 0.00445 | 2813 |
| all codes | 8 | n, k | 0.6978 | 0.00341 | 2813 |
| all codes | 8 | n, delta (theory) | 0.0529 | 0.00604 | 2813 |
| all codes | 8 | n, k + delta | 0.9553 | 0.00131 | 2813 |
| all codes | pooled | n only | 0.2157 | 0.00963 | 34775 |
| all codes | pooled | n, k | 0.4225 | 0.00827 | 34775 |
| all codes | pooled | n, delta (theory) | 0.0481 | 0.01061 | 34775 |
| all codes | pooled | n, k + delta | 0.5680 | 0.00715 | 34775 |
| all codes | pooled | n, q | 0.3714 | 0.00863 | 34775 |
| all codes | pooled | n, k, q | 0.5426 | 0.00736 | 34775 |
| all codes | pooled | n, k, q + delta | 0.7123 | 0.00584 | 34775 |
| below entropy peak | 2 | n only | 0.9062 | 0.00323 | 13404 |
| below entropy peak | 2 | n, k | 0.9535 | 0.00228 | 13404 |
| below entropy peak | 2 | n, delta (theory) | 0.9287 | 0.00282 | 13404 |
| below entropy peak | 2 | n, k + delta | 0.9568 | 0.00219 | 13404 |
| below entropy peak | 3 | n only | 0.8343 | 0.00325 | 3965 |
| below entropy peak | 3 | n, k | 0.9202 | 0.00226 | 3965 |
| below entropy peak | 3 | n, delta (theory) | 0.9189 | 0.00227 | 3965 |
| below entropy peak | 3 | n, k + delta | 0.9522 | 0.00175 | 3965 |
| below entropy peak | 4 | n only | 0.8509 | 0.00269 | 3210 |
| below entropy peak | 4 | n, k | 0.9462 | 0.00162 | 3210 |
| below entropy peak | 4 | n, delta (theory) | 0.9635 | 0.00133 | 3210 |
| below entropy peak | 4 | n, k + delta | 0.9754 | 0.00109 | 3210 |
| below entropy peak | 5 | n only | 0.8586 | 0.00256 | 2519 |
| below entropy peak | 5 | n, k | 0.9655 | 0.00127 | 2519 |
| below entropy peak | 5 | n, delta (theory) | 0.9838 | 0.00087 | 2519 |
| below entropy peak | 5 | n, k + delta | 0.9873 | 0.00077 | 2519 |
| below entropy peak | 7 | n only | 0.9089 | 0.00194 | 2729 |
| below entropy peak | 7 | n, k | 0.9837 | 0.00082 | 2729 |
| below entropy peak | 7 | n, delta (theory) | 0.9940 | 0.00050 | 2729 |
| below entropy peak | 7 | n, k + delta | 0.9982 | 0.00028 | 2729 |
| below entropy peak | 8 | n only | 0.8768 | 0.00204 | 2553 |
| below entropy peak | 8 | n, k | 0.9726 | 0.00096 | 2553 |
| below entropy peak | 8 | n, delta (theory) | 0.9934 | 0.00047 | 2553 |
| below entropy peak | 8 | n, k + delta | 0.9969 | 0.00033 | 2553 |
| below entropy peak | pooled | n only | 0.4923 | 0.00641 | 28380 |
| below entropy peak | pooled | n, k | 0.6323 | 0.00546 | 28380 |
| below entropy peak | pooled | n, delta (theory) | 0.9250 | 0.00247 | 28380 |
| below entropy peak | pooled | n, k + delta | 0.9530 | 0.00195 | 28380 |
| below entropy peak | pooled | n, q | 0.8916 | 0.00296 | 28380 |
| below entropy peak | pooled | n, k, q | 0.9484 | 0.00204 | 28380 |
| below entropy peak | pooled | n, k, q + delta | 0.9586 | 0.00183 | 28380 |
| past entropy peak | 2 | n only | 0.0088 | 0.01497 | 4478 |
| past entropy peak | 2 | n, k | 0.1424 | 0.01392 | 4478 |
| past entropy peak | 2 | n, delta (theory) | 0.0429 | 0.01471 | 4478 |
| past entropy peak | 2 | n, k + delta | 0.3887 | 0.01175 | 4478 |
| past entropy peak | 3 | n only | 0.0018 | 0.01269 | 546 |
| past entropy peak | 3 | n, k | 0.1252 | 0.01188 | 546 |
| past entropy peak | 3 | n, delta (theory) | 0.0379 | 0.01246 | 546 |
| past entropy peak | 3 | n, k + delta | 0.4256 | 0.00963 | 546 |
| past entropy peak | 4 | n only | 0.0242 | 0.00718 | 497 |
| past entropy peak | 4 | n, k | 0.0706 | 0.00700 | 497 |
| past entropy peak | 4 | n, delta (theory) | 0.0372 | 0.00713 | 497 |
| past entropy peak | 4 | n, k + delta | 0.4622 | 0.00533 | 497 |
| past entropy peak | 5 | n only | 0.0308 | 0.00433 | 372 |
| past entropy peak | 5 | n, k | 0.0674 | 0.00425 | 372 |
| past entropy peak | 5 | n, delta (theory) | 0.0319 | 0.00433 | 372 |
| past entropy peak | 5 | n, k + delta | 0.6012 | 0.00278 | 372 |
| past entropy peak | 7 | n only | 0.0005 | 0.00178 | 242 |
| past entropy peak | 7 | n, k | 0.0858 | 0.00171 | 242 |
| past entropy peak | 7 | n, delta (theory) | 0.0057 | 0.00178 | 242 |
| past entropy peak | 7 | n, k + delta | 0.7611 | 0.00087 | 242 |
| past entropy peak | 8 | n only | 0.0466 | 0.00153 | 260 |
| past entropy peak | 8 | n, k | 0.0608 | 0.00152 | 260 |
| past entropy peak | 8 | n, delta (theory) | 0.0139 | 0.00155 | 260 |
| past entropy peak | 8 | n, k + delta | 0.7011 | 0.00086 | 260 |
| past entropy peak | pooled | n only | 0.0006 | 0.01337 | 6395 |
| past entropy peak | pooled | n, k | 0.1125 | 0.01260 | 6395 |
| past entropy peak | pooled | n, delta (theory) | 0.0331 | 0.01315 | 6395 |
| past entropy peak | pooled | n, k + delta | 0.2131 | 0.01186 | 6395 |
| past entropy peak | pooled | n, q | 0.0138 | 0.01328 | 6395 |
| past entropy peak | pooled | n, k, q | 0.1220 | 0.01253 | 6395 |
| past entropy peak | pooled | n, k, q + delta | 0.3675 | 0.01064 | 6395 |
