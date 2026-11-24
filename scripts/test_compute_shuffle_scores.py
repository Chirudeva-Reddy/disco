"""Self-check of the nested checkpoints. Run: python scripts/test_compute_shuffle_scores.py path/to/dataset/"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import compute_shuffle_scores as c

tasks = c.make_tasks(Path(sys.argv[1]), 1000, None)
task = tasks[0]

nested = c.score_tf(task._replace(shuffles_per_seed=1000, seeds=c.SEEDS[:1], checkpoints=(200, 1000)))
plain = c.score_tf(task._replace(shuffles_per_seed=1000, seeds=c.SEEDS[:1]))

at_1000 = nested[nested["K"] == 1000].reset_index(drop=True)
assert at_1000.equals(plain), "K=1000 checkpoint differs from a plain 1,000-shuffle run"
at_200 = nested[nested["K"] == 200].reset_index(drop=True)
assert np.all(at_200["Shuffles_Reaching"].to_numpy() <= at_1000["Shuffles_Reaching"].to_numpy()), "reaching decreased with K"
assert set(c.score_tf(task._replace(shuffles_per_seed=100))["K"]) == {200}, "default run should report once at 2 x 100"
# Checkpoints spanning both seeds end on the default two-seed run.
both = c.make_tasks(Path(sys.argv[1]), 300, [task.genus_num], (100, 600))[0]
assert both.seeds == c.SEEDS and both.shuffles_per_seed == 300
spanning = c.score_tf(both)
assert spanning[spanning["K"] == 600].reset_index(drop=True).equals(c.score_tf(task._replace(shuffles_per_seed=300))), \
    "K=600 across two seeds differs from the default run"
print("ok")
