from math import ceil

import streamlit as st
import pandas as pd
import plotly.express as px

from app.utils import data_loading
from app.utils import constants
from app.utils import helper
from app.utils import patterns

#region Config
st.set_page_config(
    page_title=f"Shuffle K Comparison | {constants.APP_NAME}",
    initial_sidebar_state="expanded",
)

st.html(f"""
<h1 align="center" style="font-size: 48px;">
    Shuffle K Comparison |
    <img src="{constants.ASSET_LOGO_DATAURI}" alt="" width="64" height="64"/>
    {constants.APP_NAME}
</h1>
""")
st.caption("How the shuffle rank of the ELM classes within a TF changes with the number of shuffles K", text_alignment="center")

#endregion

#region Data
data = data_loading.init()
if data.shuffle_scores_by_k_df is None:
    st.warning(f"This dataset has no `{constants.PATH_DATA.TF_SHUFFLE_SCORES_BY_K.name}`. Create it with `scripts/compute_shuffle_scores.py --checkpoints 10000 50000 100000 200000`.")
    st.stop()

@st.cache_data
def compute(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ranked = patterns.rank_by_k(scores)
    return ranked, patterns.k_tf_stats(ranked), patterns.k_pairs(ranked)

ranked, tf_stats, pairs = compute(data.shuffle_scores_by_k_df)
ALL_K = sorted(ranked["K"].unique())
LARGEST_K = ALL_K[-1]
PAIRS = list(pairs[["K_A", "K_B"]].drop_duplicates().itertuples(index=False, name=None))

elm_ids = data.matches_df.drop_duplicates("ELM_Acc").set_index("ELM_Acc")["ELM_Id"]
def tf_name(genus_num: str) -> str:
    """Display name of a TF, e.g. `"P49639 | 3.1.1.2.2 | Hox-B2"`."""
    return data.genus_num_name_map.get(genus_num, genus_num)

def short_k(k: int) -> str:
    """Short label of a K, e.g. `10000` -> `"10k"`."""
    return f"{k // 1000}k"

def pair_label(k_a: int, k_b: int) -> str:
    """Label of a pair of K, e.g. `"10k vs 200k"`."""
    return f"{short_k(k_a)} vs {short_k(k_b)}"

pairs["Pair"] = [pair_label(a, b) for a, b in zip(pairs["K_A"], pairs["K_B"])]
pair_order = [pair_label(*pair) for pair in PAIRS]

#endregion

#region How to read
with st.expander("How to read this page", icon=":material/menu_book:", expanded=True):
    helper.container_indent().markdown(f"""
Every TF has its own list of ELM classes, ranked by **Shuffle p** (rank 1 =
smallest p). The ranking is worked out four times, after
{", ".join(f"{k:,}" for k in ALL_K)} shuffles ($K$). Each comparison below takes
**two** of these rankings of the **same TF** and asks how alike they are.

:primary[**Top 10%**] — the best-ranked tenth of a TF's classes. A TF with 71
classes has a top 10% of $\\lceil 7.1 \\rceil = 8$ classes. If several classes
tie at the cut, all of them are kept.

:primary[**Top-10% overlap**] — how many classes the two top-10% lists share,
out of all classes that appear in either list
([Jaccard index](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x)).
*Example:* the top 10% at 10k is {{A, B, C, D, E, F, G, H}}, at 100k it is
{{A, B, C, D, E, F, G, J}}. They share 7 classes; 9 classes appear in either
list, so the overlap is 7 / 9 = 0.78. **1.00** = exactly the same top 10%;
**0** = nothing in common.

:primary[**Spearman ρ**] — how alike the **whole** ranking is, not only the top
10% ([Spearman, 1904](https://doi.org/10.2307/1412159)). **1** = the same order
for every class, **0** = the two orders are unrelated, **−1** = reversed.
Swapping two neighbouring classes in a list of 70 leaves ρ above 0.9999,
so values here are shown with 4 decimals.

:primary[**Median**] — the middle value when the 1-per-TF numbers are sorted:
half of the TFs are at or above it. A median overlap of 1.00 means **at least
half** of the TFs keep exactly the same top 10%, so it hides the TFs that do
change. That is why the page also shows the **mean** (average), the
**lowest** value (the worst TF), and the **share of TFs with an identical top
10%**.

:primary[**At floor**] — classes that no shuffle matched as often as the real TF.
They all get the smallest possible $p = 1/(K+1)$ and therefore tie. More
shuffles can separate them, which is one reason a ranking changes with $K$.
""")

with st.expander("Method and sources", icon=":material/science:"):
    helper.container_indent().markdown(f"""
For each TF, residue 1 is kept and residues 2 onwards are shuffled
([Wregex 3.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC11318550/)). After each
shuffle, the non-overlapping matches of every class are counted. At $K$
shuffles, $p = (1 + b)/(K + 1)$, where $b$ is the number of shuffles with at
least as many matches as the real TF
([Phipson and Smyth, 2010](https://doi.org/10.2202/1544-6115.1585)).

The $K$ are **nested**: one stream of {LARGEST_K:,} shuffles per TF (two seeds of
100,000, one after the other) is read off at {", ".join(f"{k:,}" for k in ALL_K)}.
The smaller $K$ are the first shuffles of the larger, so the page shows how the
ranking settles as shuffles accumulate, not seed-to-seed variation. At 200,000
the values equal the original two-seed run.

Many shuffles per sequence to count ELM matches follow
[mimicINT](https://pmc.ncbi.nlm.nih.gov/articles/PMC11982803/) and
[Hagai et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4089993/).
Ordering classes within a TF, and comparing $K$, is our application of these
sources. A high rank says the count is unusual under reordering, not that the
match is functional.
""")

#endregion

#region All TFs
n_tfs = pairs["Genus_Num"].nunique()
st.header(f":material/insights: All {n_tfs:,} TFs")

summary = pairs.groupby("Pair", sort=False).agg(
    Same=("Top_Jaccard", lambda x: (x == 1).mean()),
    Changed=("Top_Jaccard", lambda x: int((x < 1).sum())),
    Mean_J=("Top_Jaccard", "mean"),
    Median_J=("Top_Jaccard", "median"),
    Min_J=("Top_Jaccard", "min"),
    Median_R=("Spearman", "median"),
    Min_R=("Spearman", "min"),
).reindex(pair_order).reset_index()

st.write("Each row compares two values of K, across all TFs:")
st.dataframe(
    summary,
    hide_index=True,
    column_config={
        "Pair": st.column_config.TextColumn("Comparison"),
        "Same": st.column_config.ProgressColumn("TFs with identical top 10%", min_value=0, max_value=1, format="percent"),
        "Changed": st.column_config.NumberColumn("TFs whose top 10% changed"),
        "Mean_J": st.column_config.NumberColumn("Mean top-10% overlap", format="%.3f"),
        "Median_J": st.column_config.NumberColumn("Median top-10% overlap", format="%.2f"),
        "Min_J": st.column_config.NumberColumn("Lowest top-10% overlap", format="%.2f"),
        "Median_R": st.column_config.NumberColumn("Median Spearman ρ", format="%.4f"),
        "Min_R": st.column_config.NumberColumn("Lowest Spearman ρ", format="%.4f"),
    },
)

def pair_matrix(values: pd.Series) -> pd.DataFrame:
    """Square K x K table of a per-pair value, mirrored, with 1 on the diagonal."""
    m = pd.DataFrame(1.0, index=[short_k(k) for k in ALL_K], columns=[short_k(k) for k in ALL_K])
    for (k_a, k_b), v in zip(PAIRS, values):
        m.loc[short_k(k_a), short_k(k_b)] = m.loc[short_k(k_b), short_k(k_a)] = v
    return m

left, right = st.columns(2)
for column, values, title, fmt in [
    (left, summary["Same"], "Share of TFs with identical top 10%", ".0%"),
    (right, summary["Mean_J"], "Mean top-10% overlap", ".3f"),
]:
    fig = px.imshow(pair_matrix(values), text_auto=fmt, zmin=values.min(), zmax=1,
        color_continuous_scale="Blues", labels={"x": "K", "y": "K", "color": ""}, title=title)
    fig.update_layout(coloraxis_showscale=False)
    column.plotly_chart(fig, width="stretch")

fig = px.strip(pairs, x="Pair", y="Top_Jaccard", hover_data=["Genus_Num"],
    category_orders={"Pair": pair_order},
    labels={"Pair": "", "Top_Jaccard": "Top-10% overlap"},
    title="Top-10% overlap of every TF, per comparison (one dot per TF)")
st.plotly_chart(fig, width="stretch")

wide = pairs.pivot(index="Genus_Num", columns="Pair", values=["Top_Jaccard", "Spearman"])
floor = tf_stats.pivot(index="Genus_Num", columns="K", values="At_Floor")
first = pair_label(ALL_K[0], LARGEST_K)
table = pd.DataFrame({
    "Genus_Num": wide.index,
    "TF": wide.index.map(tf_name),
    "Classes": tf_stats[tf_stats["K"] == LARGEST_K].set_index("Genus_Num").loc[wide.index, "Classes"].to_numpy(),
    **{f"J {p}": wide[("Top_Jaccard", p)].to_numpy() for p in pair_order},
    **{f"ρ {p}": wide[("Spearman", p)].to_numpy() for p in pair_order},
    **{f"Floor {short_k(k)}": floor.loc[wide.index, k].to_numpy() for k in ALL_K},
}).sort_values([f"J {first}", f"ρ {first}"], ignore_index=True)

st.write(f"One row per TF, least stable ({first}) first. Select a row to inspect it below.")
table_state = st.dataframe(
    table,
    hide_index=True,
    on_select="rerun",
    selection_mode="single-row",
    column_config={
        **{f"J {p}": st.column_config.ProgressColumn(f"Overlap {p}", min_value=0, max_value=1, format="%.2f") for p in pair_order},
        **{f"ρ {p}": st.column_config.NumberColumn(f"ρ {p}", format="%.4f") for p in pair_order},
        **{f"Floor {short_k(k)}": st.column_config.NumberColumn(f"At floor {short_k(k)}") for k in ALL_K},
    },
)
selected_rows = helper.get_df_selected_rows(table_state)

#endregion

#region Single TF
st.header(":material/visibility: One TF")

genus_options = sorted(tf_stats["Genus_Num"].unique())
selected = st.selectbox(
    "TF",
    genus_options,
    index=genus_options.index(table.loc[selected_rows[0], "Genus_Num"]) if selected_rows else 0,
    format_func=tf_name,
)

s = tf_stats[tf_stats["Genus_Num"] == selected].set_index("K")
for column, k in zip(st.columns(len(ALL_K)), ALL_K):
    column.metric(f"K = {k:,}: in top 10%", s.loc[k, "Top_Size"],
        f"{s.loc[k, 'At_Floor']} at floor p", delta_color="off", delta_arrow="off")

tf_pairs = pairs[pairs["Genus_Num"] == selected]
st.dataframe(
    tf_pairs[["Pair", "Top_Jaccard", "Spearman"]],
    hide_index=True,
    column_config={
        "Pair": st.column_config.TextColumn("Comparison"),
        "Top_Jaccard": st.column_config.ProgressColumn("Top-10% overlap", min_value=0, max_value=1, format="%.2f"),
        "Spearman": st.column_config.NumberColumn("Spearman ρ", format="%.4f"),
    },
)

tf = ranked[ranked["Genus_Num"] == selected]
tf_top = tf[tf["ELM_Acc"].isin(tf.loc[tf["Top"], "ELM_Acc"].unique())]
detail = tf_top.pivot(index="ELM_Acc", columns="K", values=["Shuffle_Rank", "Shuffle_P", "Top"])
detail_df = pd.DataFrame({
    "ELM_Acc": detail.index,
    "ELM_Id": detail.index.map(elm_ids),
    "Observed": tf_top.drop_duplicates("ELM_Acc").set_index("ELM_Acc").loc[detail.index, "Observed"].to_numpy(),
    **{f"Rank {short_k(k)}": detail[("Shuffle_Rank", k)].to_numpy() for k in ALL_K},
    **{f"p {short_k(k)}": detail[("Shuffle_P", k)].to_numpy() for k in ALL_K},
    "In top 10% at": [", ".join(short_k(k) for k in ALL_K if detail.loc[acc, ("Top", k)]) for acc in detail.index],
}).sort_values(f"Rank {short_k(LARGEST_K)}", ignore_index=True)

st.write(f"Classes in the top 10% at any K (:primary[**{len(detail_df)}**] of {s.loc[LARGEST_K, 'Classes']}):")
st.dataframe(
    detail_df,
    hide_index=True,
    column_config={
        "ELM_Acc": st.column_config.TextColumn("ELM accession"),
        "ELM_Id": st.column_config.TextColumn("ELM identifier"),
        **{f"p {short_k(k)}": st.column_config.NumberColumn(f"p {short_k(k)}", format="%.6f") for k in ALL_K},
    },
)

slope = tf_top.assign(K=lambda d: d["K"].map(lambda k: f"{k:,}"), ELM_Id=lambda d: d["ELM_Acc"].map(elm_ids))
fig = px.line(slope, x="K", y="Shuffle_Rank", color="ELM_Id", markers=True,
    category_orders={"K": [f"{k:,}" for k in ALL_K]},
    hover_data=["ELM_Acc", "Shuffle_P", "Shuffles_Reaching"],
    labels={"K": "Number of shuffles K", "Shuffle_Rank": "Shuffle rank (1 = smallest p)", "ELM_Id": "ELM class"})
fig.update_xaxes(type="category")
fig.update_yaxes(autorange="reversed")
fig.add_hline(y=ceil(patterns.TOP_FRACTION * s.loc[LARGEST_K, "Classes"]) + 0.5, line_dash="dot",
    annotation_text="top-10% cut", annotation_position="bottom right")
st.plotly_chart(fig, width="stretch")

#endregion
