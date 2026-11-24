"""
Precompute the shuffle p-value of every ELM class within every TF.

For each TF, residue 1 is kept in place and residues 2 to L are randomly
reordered. The non-overlapping regex matches of each class are then counted in
the shuffled sequence. The p-value is the share of shuffles which match at least
as often as the real sequence: `p = (1 + b) / (K + 1)`.

Sources:
- Keeping residue 1 and shuffling the rest: Wregex 3.0 (Prieto et al. 2024).
- `(1 + b) / (K + 1)`: Phipson and Smyth 2010, Hemerik and Goeman 2018.
- Many shuffles per sequence, counting ELM hits: mimicINT (2025), Hagai et al. 2014.

Using this shuffle to order the ELM classes within one TF is our application of
the above sources.

With `--checkpoints`, a single stream of shuffles (the seeds' blocks one after
another) is read off at each given K, so that the ranking at K = 10,000, 50,000,
100,000 and 200,000 can be compared (nested design: the smaller K are the first
shuffles of the larger). Its first 100,000 shuffles are the first seed's block,
and at 200,000 it is the same as the default run.
Comparing K this way is also our application.

Usage (from the repository root):
    python scripts/compute_shuffle_scores.py path/to/dataset/
    python scripts/compute_shuffle_scores.py path/to/dataset/ --genus 3.1.1.2.2 --out hoxb2.tsv
    python scripts/compute_shuffle_scores.py path/to/dataset/ --checkpoints 10000 50000 100000 200000
"""

import argparse
import math
import multiprocessing as mp
import re
import time
from pathlib import Path
from typing import Final, NamedTuple

import numpy as np
import pandas as pd


#region Constants
SEEDS: Final = (20261109, 20261110)
"""One seed per block of shuffles. The blocks are pooled."""

SHUFFLES_PER_SEED: Final = 100_000
"""mimicINT's default number of shuffles. Two seeds give K = 200,000, the value
reported for Hox-B2 in our comparison. This has not been derived for the whole
catalogue."""

OUTPUT_FILENAME: Final = "elm_shuffle_scores.tsv"
OUTPUT_FILENAME_BY_K: Final = "elm_shuffle_scores_by_k.tsv"
#endregion

#region Shuffling
class TFTask(NamedTuple):
    """The work for a single TF."""
    genus_num: str
    genus_index: int
    """Position of this TF in the sorted list of TFs, so that each TF gets its
    own reproducible random stream."""
    sequence: str
    elm_accs: list[str]
    regexes: list[str]
    observed: list[int]
    shuffles_per_seed: int
    seeds: tuple[int, ...] = SEEDS
    checkpoints: tuple[int, ...] = ()
    """Total shuffle counts (across seeds) at which to report. Empty means only
    at the end."""

def score_tf(task: TFTask) -> pd.DataFrame:
    """Compute the shuffle statistics of each class which matches the given TF.

    :param task: The TF, along with its matching classes and their real counts.
    :return: One row per class, as described in `compute_all()`.
    """

    # ! Must be the same non-overlapping scan which produced `matches.tsv`, so
    # that the real sequence and the shuffles are counted in the same way.
    compiled = [re.compile(regex.encode()) for regex in task.regexes]
    observed = np.array(task.observed, dtype=np.int64)

    real = np.frombuffer(task.sequence.encode(), dtype=np.uint8)
    real_counts = np.array([len(c.findall(real.tobytes())) for c in compiled], dtype=np.int64)
    if not np.array_equal(real_counts, observed):
        raise ValueError(f"{task.genus_num}: regex counts do not match matches.tsv")

    reaching = np.zeros(len(compiled), dtype=np.int64)
    total = np.zeros(len(compiled), dtype=np.int64)
    checkpoints = set(task.checkpoints or [len(task.seeds) * task.shuffles_per_seed])
    rows: list[pd.DataFrame] = []
    done = 0

    for seed in task.seeds:
        rng = np.random.default_rng([seed, task.genus_index])
        shuffled = real.copy()
        tail = shuffled[1:] # a view, so that shuffling it leaves residue 1 in place
        for _ in range(task.shuffles_per_seed):
            rng.shuffle(tail)
            sequence = shuffled.tobytes()
            counts = np.fromiter((len(c.findall(sequence)) for c in compiled), dtype=np.int64, count=len(compiled))
            reaching += counts >= observed
            total += counts
            done += 1
            if done in checkpoints:
                rows.append(pd.DataFrame({
                    "Genus_Num": task.genus_num,
                    "ELM_Acc": task.elm_accs,
                    "Observed": observed,
                    "Shuffle_Mean": total / done,
                    "Shuffles_Reaching": reaching.copy(),
                    "K": done,
                    "Shuffle_P": (1 + reaching) / (done + 1),
                }))

        if shuffled[0] != real[0] or not np.array_equal(np.sort(shuffled), np.sort(real)):
            raise ValueError(f"{task.genus_num}: shuffle changed residue 1 or the composition")

    if len(rows) != len(checkpoints):
        raise ValueError(f"{task.genus_num}: checkpoints beyond the number of shuffles")
    return pd.concat(rows, ignore_index=True)

#endregion

#region Dataset
def read_fasta(path: Path) -> dict[str, str]:
    """Read the FASTA file into a dictionary of Genus Number -> Sequence."""
    sequences: dict[str, list[str]] = {}
    current: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = sequences.setdefault(line[1:].partition("_")[0], [])
        else:
            current.append(line.strip())
    return {genus_num: "".join(parts) for genus_num, parts in sequences.items()}

def make_tasks(dataset_dir: Path, shuffles_per_seed: int, only_genus: list[str] | None, checkpoints: tuple[int, ...] = ()) -> list[TFTask]:
    """Build one task per TF from the dataset files.

    :param dataset_dir: The directory holding `matches.tsv` and `tfs_fasta/`.
    :param shuffles_per_seed: Number of shuffles for each seed.
    :param only_genus: If given, restrict the run to these Genus Numbers.
    :param checkpoints: If given, use as many seeds as needed to reach the
        largest checkpoint, and report at each of them.
    """
    sequences = read_fasta(dataset_dir / "tfs_fasta/tfs_all.fasta")
    matches = pd.read_csv(dataset_dir / "matches.tsv", sep="\t", usecols=["Genus_Num", "ELM_Acc", "Regex"])
    # same filter as `__load_matches()` in `app/utils/data_loading.py`
    matches = matches[matches.groupby("Regex")["ELM_Acc"].transform("count") > 1]
    per_tf = matches.groupby(["Genus_Num", "ELM_Acc", "Regex"]).size().rename("Observed").reset_index()

    n_seeds = math.ceil(max(checkpoints) / shuffles_per_seed) if checkpoints else len(SEEDS)
    if n_seeds > len(SEEDS):
        raise ValueError(f"checkpoints need {n_seeds} seeds, only {len(SEEDS)} defined")
    if checkpoints and n_seeds == 1:
        shuffles_per_seed = max(checkpoints) # stop at the last checkpoint

    tasks: list[TFTask] = []
    for genus_index, genus_num in enumerate(sorted(sequences)):
        if only_genus and genus_num not in only_genus: continue
        rows = per_tf[per_tf["Genus_Num"] == genus_num]
        if rows.empty: continue
        tasks.append(TFTask(
            genus_num=genus_num,
            genus_index=genus_index,
            sequence=sequences[genus_num],
            elm_accs=rows["ELM_Acc"].tolist(),
            regexes=rows["Regex"].tolist(),
            observed=rows["Observed"].tolist(),
            shuffles_per_seed=shuffles_per_seed,
            seeds=SEEDS[:n_seeds],
            checkpoints=checkpoints,
        ))
    return tasks

#endregion

def compute_all(tasks: list[TFTask], workers: int) -> pd.DataFrame:
    """Score all the given TFs in parallel.

    :return: One row per (TF, class) with columns: `Genus_Num`, `ELM_Acc`,
        `Observed`, `Shuffle_Mean`, `Shuffles_Reaching`, `K`, `Shuffle_P`.
    """
    started = time.perf_counter()
    results: list[pd.DataFrame] = []
    with mp.get_context("spawn").Pool(workers) as pool:
        for done, result in enumerate(pool.imap_unordered(score_tf, tasks), start=1):
            results.append(result)
            print(f"{done}/{len(tasks)} TFs, {time.perf_counter() - started:.0f} s", flush=True)
    return pd.concat(results, ignore_index=True).sort_values(["Genus_Num", "K", "ELM_Acc"], ignore_index=True)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset_dir", type=Path, help="Directory of the TF-DISCO dataset.")
    parser.add_argument("--out", type=Path, help=f"Output TSV. Default: <dataset_dir>/{OUTPUT_FILENAME}")
    parser.add_argument("--genus", nargs="*", help="Only score these Genus Numbers.")
    parser.add_argument("--shuffles-per-seed", type=int, default=SHUFFLES_PER_SEED)
    parser.add_argument("--checkpoints", type=int, nargs="*", default=[],
        help=f"Report at each of these K, from one stream (first seed). Default output: <dataset_dir>/{OUTPUT_FILENAME_BY_K}")
    parser.add_argument("--workers", type=int, default=mp.cpu_count())
    args = parser.parse_args()

    checkpoints = tuple(sorted(args.checkpoints))
    tasks = make_tasks(args.dataset_dir, args.shuffles_per_seed, args.genus, checkpoints)
    scores = compute_all(tasks, args.workers)
    out: Path = args.out or (args.dataset_dir / (OUTPUT_FILENAME_BY_K if checkpoints else OUTPUT_FILENAME))
    scores.to_csv(out, sep="\t", index=False)
    print(f"Wrote {len(scores)} rows to {out}", flush=True)

if __name__ == "__main__":
    main()
