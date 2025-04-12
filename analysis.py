import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

# Updated optimization data.
optimizations = OrderedDict([
    ("Baseline", {
        "1000/5000": [27.154503, 27.282169, 27.287363],
        "5/1000000000": [168.969376, 171.720993, 168.426620],
        "1000/100000": [542.263000, 540.893071, 543.775002]
    }),
    ("Clang++", {
        "1000/5000": [25.629236, 25.581055, 25.583626],
        "5/1000000000": [162.166168, 161.630447, 161.376495],
        "1000/100000": [511.597870, 510.129103, 511.912018]
    }),
    ("Memcpy, Mass caching, exploiting symmetry", {
        "1000/5000": [14.821008, 14.756897, 14.792585],
        "5/1000000000": [118.436707, 121.182762, 118.671623],
        "1000/100000": [296.270081, 298.239821, 296.005399]
    }),
    ("inline next, compiler flags", {
        "1000/5000": [12.665927, 12.434470, 12.556138],
        "5/1000000000": [79.620979, 77.496124, 77.580162],
        "1000/100000": [251.394531, 250.097184, 251.577190]
    }),
    ("Reducing mallocs/frees", {
        "1000/5000": [12.209402, 12.230252, 12.250868],
        "5/1000000000": [67.463959, 69.044136, 69.372330],
        "1000/100000": [242.839752, 240.309101, 243.209589]
    }),
    ("Parallelize", {
        "1000/5000": [5.354763, 5.155138, 5.621915],
        "5/1000000000": [70.149384, 69.751678, 69.791748],
        "1000/100000": [107.741165, 112.062477, 109.468536]
    })
])

# Update datasets to reflect the new keys.
datasets = ["1000/5000", "5/1000000000", "1000/100000"]

ordered_labels = list(optimizations.keys())
n_steps = len(ordered_labels)
transition_labels = [ordered_labels[i] for i in range(1, n_steps)]

# Calculate relative speedups (previous step’s min divided by current step’s min)
speedups = {ds: [] for ds in datasets}
for i in range(1, n_steps):
    prev_opt = optimizations[ordered_labels[i - 1]]
    curr_opt = optimizations[ordered_labels[i]]
    for ds in datasets:
        prev_vals = prev_opt.get(ds)
        curr_vals = curr_opt.get(ds)
        prev_min = np.min(prev_vals)
        curr_min = np.min(curr_vals)
        speedups[ds].append(prev_min / curr_min)

for ds in datasets:
    speedups[ds] = np.array(speedups[ds])

# Plotting
x = np.arange(len(transition_labels))
width = 0.25

fig, ax = plt.subplots(figsize=(14, 6))

# Define offsets and colors for each dataset.
offsets = {
    "1000/5000": -width,
    "5/1000000000": 0,
    "1000/100000": width
}
colors = {
    "1000/5000": 'skyblue',
    "5/1000000000": 'lightgreen',
    "1000/100000": 'salmon'
}

for ds in datasets:
    ax.bar(x + offsets[ds], speedups[ds], width, label=ds, color=colors[ds])

ax.set_xlabel("Optimization Transition")
ax.set_ylabel("Relative Speedup (x)")
ax.set_xticks(x)
ax.set_xticklabels(transition_labels, rotation=45, ha="right")
ax.legend()

# Annotate bars with speedup values.
for ds in datasets:
    for i, val in enumerate(speedups[ds]):
        ax.text(x[i] + offsets[ds], val + 0.02, f"{val:.2f}x", 
                ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.show()
plt.savefig('optimization_results.png', dpi=600, bbox_inches='tight')
