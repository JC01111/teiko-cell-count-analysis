"""Interactive dashboard for the Loblaw Bio cell-count dataset.

Run with:  streamlit run dashboard/app.py     (or `make dashboard`)

Everything the dashboard shows is read from the SQLite database built by
`load_data.py`; if the database is missing it is built on first launch so the
app also works on a fresh clone / a cloud deployment.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from analysis.common import DB_PATH, POPULATION_LABELS, analysis_frame, population_order  # noqa: E402
from analysis.statistics_report import run_tests  # noqa: E402

RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}
COLORS = {"Responder": "#2E7D9A", "Non-responder": "#C05746"}

st.set_page_config(page_title="Loblaw Bio | Cell Count Explorer", page_icon="🧬",
                   layout="wide", initial_sidebar_state="expanded")


# --------------------------------------------------------------------- data
@st.cache_data(show_spinner="Loading data…")
def get_data() -> pd.DataFrame:
    if not os.path.exists(DB_PATH):
        import load_data
        load_data.main([])
    df = analysis_frame()
    df["response_label"] = df["response"].map(RESPONSE_LABELS).fillna("Unknown")
    return df


def multiselect_all(label: str, options, key: str, default=None):
    """A multiselect that treats 'nothing chosen' as 'everything chosen'."""
    chosen = st.sidebar.multiselect(label, options, default=default if default is not None else list(options), key=key)
    return list(options) if not chosen else chosen


df = get_data()

# ------------------------------------------------------------------ sidebar
st.sidebar.title("🧬 Filters")
st.sidebar.caption("Applies to the Overview, Frequencies and Responder tabs.")

conditions = multiselect_all("Condition", sorted(df["condition"].unique()), "f_cond")
treatments = multiselect_all("Treatment", sorted(df["treatment"].unique()), "f_treat")
sample_types = multiselect_all("Sample type", sorted(df["sample_type"].unique()), "f_stype")
projects = multiselect_all("Project", sorted(df["project"].unique()), "f_proj")
sexes = multiselect_all("Sex", sorted(df["sex"].dropna().unique()), "f_sex")
times = st.sidebar.multiselect(
    "Timepoint (days from treatment start)",
    sorted(df["time_from_treatment_start"].dropna().unique().tolist()),
    default=sorted(df["time_from_treatment_start"].dropna().unique().tolist()),
    key="f_time",
)
if not times:
    times = sorted(df["time_from_treatment_start"].dropna().unique().tolist())

view = df[
    df["condition"].isin(conditions)
    & df["treatment"].isin(treatments)
    & df["sample_type"].isin(sample_types)
    & df["project"].isin(projects)
    & df["sex"].isin(sexes)
    & df["time_from_treatment_start"].isin(times)
]

st.sidebar.markdown("---")
if st.sidebar.button("Reset to the Part 3 cohort", width="stretch"):
    for k, v in {
        "f_cond": ["melanoma"], "f_treat": ["miraclib"], "f_stype": ["PBMC"],
        "f_proj": sorted(df["project"].unique()), "f_sex": sorted(df["sex"].dropna().unique()),
        "f_time": sorted(df["time_from_treatment_start"].dropna().unique().tolist()),
    }.items():
        st.session_state[k] = v
    st.rerun()

st.sidebar.caption("Database: `%s`" % os.path.relpath(DB_PATH, ROOT))

# -------------------------------------------------------------------- header
st.title("Cell Count Explorer")
st.caption(
    "Immune cell populations across clinical samples — Loblaw Bio / Teiko technical assignment."
)

k1, k2, k3, k4, k5 = st.columns(5)
k1.metric("Samples", "{:,}".format(view["sample"].nunique()))
k2.metric("Subjects", "{:,}".format(view["subject"].nunique()))
k3.metric("Projects", view["project"].nunique())
k4.metric("Populations", view["population"].nunique())
k5.metric("Total cells", "{:,.0f}".format(view.drop_duplicates("sample")["total_count"].sum()))

if view.empty:
    st.warning("No samples match the current filters.")
    st.stop()

tab_overview, tab_freq, tab_resp, tab_baseline, tab_schema = st.tabs(
    ["Overview", "Part 2 · Frequencies", "Part 3 · Responders", "Part 4 · Baseline subset", "Schema"]
)

# ------------------------------------------------------------------ overview
with tab_overview:
    left, right = st.columns([3, 2])

    with left:
        st.subheader("Population composition")
        comp = (
            view.groupby(["population_label"], as_index=False)["percentage"].mean()
            .sort_values("percentage", ascending=False)
        )
        fig = px.bar(comp, x="percentage", y="population_label", orientation="h",
                     text=comp["percentage"].map("{:.1f}%".format),
                     labels={"percentage": "Mean relative frequency (%)", "population_label": ""},
                     color_discrete_sequence=["#2E7D9A"])
        fig.update_layout(height=340, margin=dict(l=0, r=0, t=10, b=0), yaxis=dict(autorange="reversed"))
        st.plotly_chart(fig)

        st.subheader("Distribution of relative frequency")
        fig2 = px.violin(view, x="population_label", y="percentage", box=True, points=False,
                         color="population_label", labels={"population_label": "",
                                                           "percentage": "Relative frequency (%)"})
        fig2.update_layout(height=380, showlegend=False, margin=dict(l=0, r=0, t=10, b=0))
        st.plotly_chart(fig2)

    with right:
        st.subheader("Cohort make-up")
        samples = view.drop_duplicates("sample")
        for col, title in [("condition", "Condition"), ("treatment", "Treatment"),
                           ("sample_type", "Sample type"), ("response_label", "Response")]:
            counts = samples[col].value_counts().reset_index()
            counts.columns = [col, "samples"]
            fig = px.pie(counts, names=col, values="samples", hole=0.55, title=title)
            fig.update_layout(height=230, margin=dict(l=0, r=0, t=40, b=0), showlegend=True)
            fig.update_traces(textinfo="percent")
            st.plotly_chart(fig)

# --------------------------------------------------------------- frequencies
with tab_freq:
    st.subheader("Relative frequency of each cell population, per sample")
    st.caption(
        "`percentage = 100 × count / total_count`, where `total_count` is the sum of all "
        "populations measured on that sample. Computed by the `v_sample_frequency` SQL view."
    )
    table = (
        view[["sample", "total_count", "population", "count", "percentage"]]
        .sort_values(["sample", "population"])
        .reset_index(drop=True)
    )
    st.dataframe(table, width="stretch", height=430, hide_index=True)
    st.download_button(
        "Download this table as CSV",
        table.to_csv(index=False).encode("utf-8"),
        file_name="cell_frequencies.csv",
        mime="text/csv",
    )

    st.subheader("Composition of individual samples")
    n_show = st.slider("Samples to display", 5, 60, 25, step=5)
    picks = sorted(view["sample"].unique())[:n_show]
    stacked = view[view["sample"].isin(picks)]
    fig = px.bar(stacked, x="sample", y="percentage", color="population_label",
                 labels={"percentage": "Relative frequency (%)", "sample": "",
                         "population_label": "Population"})
    fig.update_layout(height=420, barmode="stack", margin=dict(l=0, r=0, t=10, b=0))
    st.plotly_chart(fig)

# ---------------------------------------------------------------- responders
with tab_resp:
    st.subheader("Responders vs non-responders")
    st.caption(
        "The assignment asks about melanoma / miraclib / PBMC — use the sidebar to point the "
        "same analysis at any other cohort."
    )
    cohort = view[view["response"].isin(["yes", "no"])].copy()
    if cohort.empty:
        st.info("The current filters contain no labelled responders or non-responders.")
    else:
        c1, c2 = st.columns([1, 1])
        c1.metric("Responder samples", "{:,}".format(
            cohort[cohort["response"] == "yes"]["sample"].nunique()))
        c2.metric("Non-responder samples", "{:,}".format(
            cohort[cohort["response"] == "no"]["sample"].nunique()))

        order = [POPULATION_LABELS.get(p, p) for p in population_order(cohort)]
        fig = px.box(cohort, x="population_label", y="percentage", color="response_label",
                     category_orders={"population_label": order,
                                      "response_label": ["Responder", "Non-responder"]},
                     color_discrete_map=COLORS, points=False,
                     labels={"population_label": "", "percentage": "Relative frequency (%)",
                             "response_label": "Response"})
        fig.update_layout(height=480, boxmode="group", margin=dict(l=0, r=0, t=10, b=0))
        st.plotly_chart(fig)

        st.markdown("**Statistical comparison** — two-sided Mann-Whitney U, "
                    "Benjamini-Hochberg FDR across populations.")
        subject_level = st.checkbox(
            "Collapse to one value per subject (removes repeated-measures pseudoreplication)",
            value=False,
        )
        data = cohort
        if subject_level:
            data = (cohort.groupby(["subject", "response", "population"], as_index=False)
                    ["percentage"].mean())
        res = run_tests(data)
        if res.empty:
            st.info("Not enough samples in one of the groups to run the tests.")
        else:
            show = res[["population_label", "n_responder", "n_non_responder",
                        "median_responder", "median_non_responder", "median_difference",
                        "rank_biserial", "p_mannwhitney", "p_adj_bh",
                        "significant_fdr_0.05", "higher_in"]]
            st.dataframe(
                show.style.format({
                    "median_responder": "{:.2f}", "median_non_responder": "{:.2f}",
                    "median_difference": "{:+.2f}", "rank_biserial": "{:+.3f}",
                    "p_mannwhitney": "{:.3g}", "p_adj_bh": "{:.3g}",
                }),
                width="stretch", hide_index=True,
            )
            hits = res[res["significant_fdr_0.05"]]
            if hits.empty:
                st.info("No population differs significantly (FDR < 0.05) in this cohort.")
            else:
                for _, r in hits.iterrows():
                    st.success(
                        "**%s** — BH p = %.3g; median %.2f%% in responders vs %.2f%% in "
                        "non-responders (%+.2f pp, higher in %s)."
                        % (r["population_label"], r["p_adj_bh"], r["median_responder"],
                           r["median_non_responder"], r["median_difference"], r["higher_in"])
                    )

        st.markdown("**Time course** — median relative frequency by days on treatment.")
        tc = (cohort.groupby(["time_from_treatment_start", "population_label", "response_label"],
                             as_index=False)["percentage"].median())
        fig = px.line(tc, x="time_from_treatment_start", y="percentage",
                      color="response_label", facet_col="population_label",
                      facet_col_wrap=5, markers=True, color_discrete_map=COLORS,
                      category_orders={"population_label": order},
                      labels={"time_from_treatment_start": "Days from treatment start",
                              "percentage": "Median frequency (%)", "response_label": "Response"})
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig.update_yaxes(matches=None, showticklabels=True)
        fig.update_layout(height=380, margin=dict(l=0, r=0, t=40, b=0))
        st.plotly_chart(fig)

# ------------------------------------------------------------------ baseline
with tab_baseline:
    st.subheader("Baseline subset")
    st.caption("Defaults to the Part 4 definition: melanoma · miraclib · PBMC · time = 0. "
               "These controls are independent of the sidebar.")

    b1, b2, b3, b4 = st.columns(4)
    b_cond = b1.selectbox("Condition", sorted(df["condition"].unique()),
                          index=sorted(df["condition"].unique()).index("melanoma"))
    b_treat = b2.selectbox("Treatment", sorted(df["treatment"].unique()),
                           index=sorted(df["treatment"].unique()).index("miraclib"))
    b_type = b3.selectbox("Sample type", sorted(df["sample_type"].unique()),
                          index=sorted(df["sample_type"].unique()).index("PBMC"))
    b_time = b4.selectbox("Timepoint", sorted(df["time_from_treatment_start"].dropna().unique()),
                          index=0)

    sub = df[(df["condition"] == b_cond) & (df["treatment"] == b_treat)
             & (df["sample_type"] == b_type) & (df["time_from_treatment_start"] == b_time)]
    samples = sub.drop_duplicates("sample")

    st.metric("Samples in subset", "{:,}".format(len(samples)))
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("**(a) Samples per project**")
        t = samples["project"].value_counts().rename_axis("project").reset_index(name="samples")
        st.dataframe(t, hide_index=True, width="stretch")
        st.plotly_chart(px.bar(t, x="project", y="samples",
                               color_discrete_sequence=["#2E7D9A"]).update_layout(
            height=260, margin=dict(l=0, r=0, t=10, b=0)))
    with c2:
        st.markdown("**(b) Responders vs non-responders**")
        t = (samples["response"].fillna("unknown").value_counts()
             .rename_axis("response").reset_index(name="samples"))
        st.dataframe(t, hide_index=True, width="stretch")
        st.plotly_chart(px.bar(t, x="response", y="samples", color="response",
                               color_discrete_sequence=["#2E7D9A", "#C05746"]).update_layout(
            height=260, showlegend=False, margin=dict(l=0, r=0, t=10, b=0)))
    with c3:
        st.markdown("**(c) Males vs females**")
        t = (samples["sex"].fillna("unknown").value_counts()
             .rename_axis("sex").reset_index(name="samples"))
        st.dataframe(t, hide_index=True, width="stretch")
        st.plotly_chart(px.bar(t, x="sex", y="samples", color="sex",
                               color_discrete_sequence=["#6C8EA4", "#C7956D"]).update_layout(
            height=260, showlegend=False, margin=dict(l=0, r=0, t=10, b=0)))

    st.markdown("---")
    st.markdown("**Average B cells — %s males, responders, time = %d, all sample types "
                "and treatments**" % (b_cond, b_time))
    bcell = df[(df["condition"] == b_cond) & (df["sex"] == "M") & (df["response"] == "yes")
               & (df["time_from_treatment_start"] == b_time) & (df["population"] == "b_cell")]
    if bcell.empty:
        st.info("No matching samples.")
    else:
        m1, m2, m3 = st.columns(3)
        m1.metric("Samples", "{:,}".format(len(bcell)))
        m2.metric("Average B cell count", "{:,.1f}".format(bcell["count"].mean()))
        m3.metric("Median B cell count", "{:,.0f}".format(bcell["count"].median()))
        fig = px.histogram(bcell, x="count", nbins=40, color_discrete_sequence=["#2E7D9A"],
                           labels={"count": "B cell count"})
        fig.add_vline(x=bcell["count"].mean(), line_dash="dash", line_color="#C05746",
                      annotation_text="mean")
        fig.update_layout(height=320, margin=dict(l=0, r=0, t=10, b=0))
        st.plotly_chart(fig)

    st.dataframe(
        samples[["sample", "subject", "project", "response", "sex", "age"]]
        .sort_values("sample").reset_index(drop=True),
        width="stretch", height=300, hide_index=True,
    )

# -------------------------------------------------------------------- schema
with tab_schema:
    st.subheader("Database schema")
    st.markdown(
        """
The CSV is normalised into a small star-ish schema. The important choice is that
**cell counts are rows, not columns** (`sample_cell_count`), so a new marker is an
`INSERT` into `cell_population` rather than a schema migration.

```
project ─┐
         ├─< enrollment >── subject
condition┤        │
treatment┘        │
                  └──< sample >──< sample_cell_count >── cell_population
                          │
                     sample_type
```

* **subject** — the person (`sex` only; nothing that changes per study).
* **enrollment** — one subject in one project on one treatment arm, carrying
  `condition`, `response` and `age_at_enrollment`. A subject can therefore appear in
  several trials without duplicating rows.
* **sample** — a specimen drawn from an enrollment at a given `time_from_treatment_start`.
* **sample_cell_count** — the fact table, one row per (sample, population).
* **v_sample_frequency / v_analysis_base** — views that define relative frequency once
  so SQL, Python and this dashboard cannot drift apart.
"""
    )
    st.code(open(os.path.join(ROOT, "db", "schema.sql")).read(), language="sql")
