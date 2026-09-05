import pandas as pd
import numpy as np

from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from sklearn.neighbors import KNeighborsRegressor

# ============================================================
# INPUT
# ============================================================

q = 2       # <-- CHANGE THIS: 2, 3, 4, 5, 7, or 8

# ============================================================
# Load combined data
# ============================================================

df = pd.read_csv(r"C:\Users\utsav\OneDrive\ドキュメント\Code Tables Data - Final\Excel Files\all_codes_combined.csv")

# ============================================================
# Select the requested q
# ============================================================

df_q = df[df["q"] == q].copy()

# ============================================================
# ONLY semi-perfect quantum codes
# ============================================================

sp = df_q[
    df_q["Type of Code"].astype(str).str.strip().str.lower()
    == "semi-perfect quantum code"
].copy()

# ============================================================
# Predictors and target
# ============================================================

X = sp[["n", "k"]]
y = sp["C'"]

print(f"q = {q}")
print("Number of semi-perfect codes:", len(sp))

# ============================================================
# Models
# ============================================================

models = {

    "KNN": make_pipeline(
        StandardScaler(),
        KNeighborsRegressor(
            n_neighbors=15,
            weights="distance"
        )
    ),

    "Random Forest": RandomForestRegressor(
        n_estimators=500,
        random_state=42,
        min_samples_leaf=5
    ),

    "Extra Trees": ExtraTreesRegressor(
        n_estimators=500,
        random_state=42,
        min_samples_leaf=5
    ),

    "Gaussian Process": make_pipeline(
        StandardScaler(),
        GaussianProcessRegressor(
            kernel=
                ConstantKernel(1.0)
                * Matern(
                    length_scale=1.0,
                    nu=1.5
                )
                + WhiteKernel(
                    noise_level=0.01
                ),
            normalize_y=True,
            random_state=42
        )
    )
}

# ============================================================
# 10-fold Cross-validation
# ============================================================

cv = KFold(
    n_splits=10,
    shuffle=True,
    random_state=42
)

results = {}

for name, model in models.items():

    errors = []

    for train_idx, test_idx in cv.split(X):

        model.fit(
            X.iloc[train_idx],
            y.iloc[train_idx]
        )

        pred = model.predict(
            X.iloc[test_idx]
        )

        errors.extend(
            np.abs(
                pred - y.iloc[test_idx]
            )
        )

    results[name] = np.mean(errors)

# ============================================================
# Results
# ============================================================

print("\n10-fold CV mean absolute error:\n")

for name, error in sorted(
    results.items(),
    key=lambda x: x[1]
):
    print(f"{name:20s}: {error:.6f}")
    
# %%

from sklearn.model_selection import GroupKFold
from sklearn.base import clone

# Group all semi-perfect codes having the same (n,k)
groups = sp[["n", "k"]].astype(str).agg("_".join, axis=1)

gcv = GroupKFold(n_splits=5)

grouped_results = {}

for name, model in models.items():

    errors = []

    for train_idx, test_idx in gcv.split(X, y, groups=groups):

        m = clone(model)

        m.fit(
            X.iloc[train_idx],
            y.iloc[train_idx]
        )

        pred = m.predict(
            X.iloc[test_idx]
        )

        errors.extend(
            np.abs(
                pred - y.iloc[test_idx]
            )
        )

    grouped_results[name] = np.mean(errors)

print(f"Grouped 5-fold CV mean absolute error (q={q}):\n")

for name, error in sorted(
    grouped_results.items(),
    key=lambda x: x[1]
):
    print(f"{name:20s}: {error:.6f}")
    
# %%
    
from scipy.special import gammaln, logsumexp
import numpy as np
from scipy.optimize import brentq
import ast
from collections import defaultdict


# ============================================================
# INPUTS
# ============================================================

n = 73
k = 18
q = 2

# Enter the w_i values here
w_values = [40,42]   # <-- CHANGE THIS

# Empirical error of C' from grouped cross-validation
C_prime_error = grouped_results[
    min(grouped_results, key=grouped_results.get)
]

# ============================================================
# 1. Calculate the summation
# ============================================================

# For general q:
# sum = Σ C(n,w_i) * (q^2 - 1)^(w_i)

log_terms = np.array([
    gammaln(n + 1)
    - gammaln(w + 1)
    - gammaln(n - w + 1)
    + w * np.log(q**2 - 1)
    for w in w_values
])

log_summation = logsumexp(log_terms)

logq_summation = log_summation / np.log(q)

# ============================================================
# 2. Calculate R
# ============================================================

R = logq_summation / (2 * n)

# ============================================================
# 3. Calculate log_10(2^(2n(n-k)) / summation^(n-k))
# ============================================================

log_ratio = (
    2 * n * (n - k) * np.log10(q)
    - (n - k) * log_summation / np.log(10)
)

# ============================================================
# 4. Get C' from the best model
# ============================================================

best_model_name = min(
    grouped_results,
    key=grouped_results.get
)

best_model = models[best_model_name]

best_model.fit(X, y)

target = pd.DataFrame({
    "n": [n],
    "k": [k]
})

C_prime = best_model.predict(target)[0]

# Account for empirical C' error
C_prime_low = C_prime - C_prime_error
C_prime_high = C_prime + C_prime_error

# ============================================================
# 5. Define the q^2-ary entropy function
# ============================================================

def Hq2(x, q):

    if x == 0:
        return 0.0

    if x == 1:
        return 1.0

    Q = q**2

    return (
        x * np.log(Q - 1) / np.log(Q)
        - x * np.log(x) / np.log(Q)
        - (1 - x) * np.log(1 - x) / np.log(Q)
    )

# ============================================================
# 6. Solve H_(q^2)(a/n) = R + C'
# ============================================================

def solve_a(C):

    target_entropy = R + C

    if target_entropy < 0 or target_entropy > 1:
        return None, target_entropy, None

    x = brentq(
        lambda x: Hq2(x, q) - target_entropy,
        0,
        (q**2 - 1) / (q**2)
    )

    a = n * x

    return a, target_entropy, x


# Central estimate
a_central, entropy_central, x_central = solve_a(C_prime)

# Lower C' -> lower a
a_low, entropy_low, x_low = solve_a(C_prime_low)

# Higher C' -> higher a
a_high, entropy_high, x_high = solve_a(C_prime_high)

# ============================================================
# 7. Weight spectrum + average weight estimator
# ============================================================

def estimate_weight_spectrum(
    n_target,
    k_target,
    num_neighbors=15
):

    # ALL codes for weight-spectrum estimation
    data = df_q.copy()

    # Distance in (n,k)
    data["distance"] = np.sqrt(
        (data["n"] - n_target)**2 +
        (data["k"] - k_target)**2
    )

    # Nearest codes
    nearest = data.sort_values("distance").head(
        num_neighbors
    )

    frequency_estimate = defaultdict(float)

    for _, row in nearest.iterrows():

        ws = row["Weight Spectrum"]

        if isinstance(ws, str):
            ws = ast.literal_eval(ws)

        distance = row["distance"]

        weight = 1.0 / (distance + 1e-6)

        for w, f in ws.items():

            frequency_estimate[int(w)] += (
                weight * int(f)
            )

    # --------------------------------------------------------
    # Convert estimated frequencies into integers
    # while forcing sum(f_i) = n-k
    # --------------------------------------------------------

    target_total = n_target - k_target

    if target_total <= 0:
        return {}

    total_estimate = sum(
        frequency_estimate.values()
    )

    if total_estimate == 0:
        return {}

    proportions = {
        w: f / total_estimate
        for w, f in frequency_estimate.items()
    }

    raw_counts = {
        w: proportions[w] * target_total
        for w in proportions
    }

    integer_counts = {
        w: int(np.floor(value))
        for w, value in raw_counts.items()
    }

    # Distribute remaining counts using largest remainder
    remaining = (
        target_total
        - sum(integer_counts.values())
    )

    remainders = sorted(
        raw_counts.keys(),
        key=lambda w:
            raw_counts[w] - integer_counts[w],
        reverse=True
    )

    for w in remainders[:remaining]:
        integer_counts[w] += 1

    # Remove zero-frequency entries
    integer_counts = {
        w: f
        for w, f in sorted(integer_counts.items())
        if f > 0
    }

    return integer_counts

# ============================================================
# 8. Output
# ============================================================

print(f"n = {n}")
print(f"k = {k}")
print(f"q = {q}")
print(f"w_i values = {w_values}")
print()

print(f"R = {R:.10f}")
print()

print(f"Best model = {best_model_name}")
print(f"C' = {C_prime:.10f}")
print(f"C' error = ±{C_prime_error:.6f}")
print(
    f"C' range = "
    f"[{C_prime_low:.10f}, {C_prime_high:.10f}]"
)
print()

print(f"R + C' = {entropy_central:.10f}")
print(
    f"R + C' range = "
    f"[{entropy_low:.10f}, {entropy_high:.10f}]"
)
print()

print(
    f"Speed Up = {log_ratio:.10f}"
)

if a_central is not None:

    print()
    print(f"a/n (central) = {x_central:.10f}")
    print(f"a (central)   = {a_central:.10f}")

    print()
    print(
        f"a/n (range) = "
        f"[{x_low:.10f}, {x_high:.10f}]"
    )

    print(
        f"a (range)   = "
        f"[{a_low:.10f}, {a_high:.10f}]"
    )

    print()
    print(
        f"a = {a_central:.10f} "
        f"+/- {(a_high - a_low) / 2:.10f}"
    )

else:

    print()
    print(
        "No solution for "
        "H_(q^2)(a/n) = R + C'."
    )

# ============================================================
# Estimated weight spectrum
# ============================================================

estimated_spectrum = estimate_weight_spectrum(
    n,
    k,
    num_neighbors=15
)

print()
print("Estimated Weight Spectrum:")
print(estimated_spectrum)

# Average weight
if estimated_spectrum:

    total_weights = sum(
        w * f
        for w, f in estimated_spectrum.items()
    )

    average_weight = (
        total_weights / (n - k)
    )

    print(
        f"Estimated Average Weight = "
        f"{average_weight:.6f}"
    )
    
#%%

# ============================================================
# Blind verification of a known (q,n,k)
# ============================================================

n_test = 50
k_test = 42
q_test = 2

# ------------------------------------------------------------
# Select semi-perfect quantum codes for this q
# ------------------------------------------------------------

sp_test = df[
    (df["q"] == q_test) &
    (
        df["Type of Code"].astype(str).str.strip().str.lower()
        == "semi-perfect quantum code"
    )
].copy()

# ------------------------------------------------------------
# Find the known semi-perfect code(s)
# ------------------------------------------------------------

target_rows = sp_test[
    (sp_test["n"] == n_test) &
    (sp_test["k"] == k_test)
]

if target_rows.empty:

    print(
        f"No semi-perfect code found for "
        f"(q,n,k)=({q_test},{n_test},{k_test})"
    )

else:

    # Actual C' value
    actual_values = target_rows["C'"].values

    # --------------------------------------------------------
    # REMOVE ALL codes with this (q,n,k) from training
    # --------------------------------------------------------

    train_mask = ~(
        (sp_test["n"] == n_test) &
        (sp_test["k"] == k_test)
    )

    X_train = sp_test.loc[
        train_mask,
        ["n", "k"]
    ]

    y_train = sp_test.loc[
        train_mask,
        "C'"
    ]

    # --------------------------------------------------------
    # Train KNN WITHOUT the target (n,k)
    # --------------------------------------------------------

    blind_model = models["KNN"]

    blind_model.fit(
        X_train,
        y_train
    )

    # --------------------------------------------------------
    # Predict the hidden value
    # --------------------------------------------------------

    target = pd.DataFrame({
        "n": [n_test],
        "k": [k_test]
    })

    predicted = blind_model.predict(target)[0]

    # --------------------------------------------------------
    # Compare with actual value
    # --------------------------------------------------------

    print(
        f"Target code: "
        f"(q,n,k)=({q_test},{n_test},{k_test})"
    )
    print()

    print(
        f"Empirical prediction = "
        f"{predicted:.10f}"
    )

    print(
        f"Actual C'            = "
        f"{actual_values[0]:.10f}"
    )

    print(
        f"Absolute error        = "
        f"{abs(predicted - actual_values[0]):.10f}"
    )