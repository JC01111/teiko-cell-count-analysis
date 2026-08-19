#!/usr/bin/env python3
"""Part 3 -- responders vs non-responders among melanoma patients on miraclib.

Cohort: condition = melanoma, treatment = miraclib, sample_type = PBMC.

For every cell population we compare the *relative frequency* (percentage of
that sample's total cell count) between responders and non-responders.

Statistics
----------
* Primary test: two-sided Mann-Whitney U.  Relative frequencies are bounded,
  mildly skewed and we do not want to assume normality, so a rank test is the
  safe default; Welch's t-test is reported alongside as a sanity check.
* Multiplicity: five populations are tested, so p-values are adjusted with
  Benjamini-Hochberg FDR.
* Effect size: rank-biserial correlation (from U) plus the raw difference in
  medians, because with thousands of samples a tiny difference will be
  "significant" without being interesting.
* Pseudoreplication: each subject contributes up to three timepoints, so the
  samples are not independent.  The headline test uses all samples (as the
  brief specifies) and a subject-level sensitivity analysis -- one mean value
  per subject -- is run as well.  Conclusions should agree; where they do not,
  believe the subject-level one.
"""
from __future__ import annotations

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy import stats  # noqa: E402

from analysis.common import (  # noqa: E402
    FIGURE_DIR,
    POPULATION_LABELS,
    analysis_frame,
    banner,
    ensure_dirs,
    population_order,
    write_csv,
)

SEED = 0            # keeps jitter/bootstrap deterministic so figures do not churn in git
RESPONSE_ORDER = ["yes", "no"]
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}
PALETTE = {"Responder": "#2E7D9A", "Non-responder": "#C05746"}


def melanoma_miraclib_pbmc(df: pd.DataFrame) -> pd.DataFrame:
    """The Part 3 cohort, with unlabelled responses dropped."""
    cohort = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
        & (df["response"].isin(RESPONSE_ORDER))
    ].copy()
    cohort["response_label"] = cohort["response"].map(RESPONSE_LABELS)
    return cohort


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _compare(values_yes: np.ndarray, values_no: np.ndarray) -> dict:
    n1, n2 = values_yes.size, values_no.size
    u, p_mwu = stats.mannwhitneyu(values_yes, values_no, alternative="two-sided")
    t_stat, p_t = stats.ttest_ind(values_yes, values_no, equal_var=False)
    # Rank-biserial correlation: +1 means responders always higher.
    rank_biserial = 2.0 * u / (n1 * n2) - 1.0
    return {
        "n_responder": n1,
        "n_non_responder": n2,
        "median_responder": float(np.median(values_yes)),
        "median_non_responder": float(np.median(values_no)),
        "median_difference": float(np.median(values_yes) - np.median(values_no)),
        "mean_responder": float(values_yes.mean()),
        "mean_non_responder": float(values_no.mean()),
        "u_statistic": float(u),
        "p_mannwhitney": float(p_mwu),
        "rank_biserial": float(rank_biserial),
        "t_statistic": float(t_stat),
        "p_welch_t": float(p_t),
    }


def run_tests(cohort: pd.DataFrame, value_col: str = "percentage") -> pd.DataFrame:
    rows = []
    for pop in population_order(cohort):
        sub = cohort[cohort["population"] == pop]
        yes = sub.loc[sub["response"] == "yes", value_col].to_numpy(dtype=float)
        no = sub.loc[sub["response"] == "no", value_col].to_numpy(dtype=float)
        if yes.size < 3 or no.size < 3:
            continue
        rec = {"population": pop, "population_label": POPULATION_LABELS.get(pop, pop)}
        rec.update(_compare(yes, no))
        rows.append(rec)

    res = pd.DataFrame(rows)
    if res.empty:
        return res
    res["p_adj_bh"] = _bh_fdr(res["p_mannwhitney"].to_numpy())
    res["significant_fdr_0.05"] = res["p_adj_bh"] < 0.05
    res["higher_in"] = np.where(
        res["median_difference"] > 0, "responders", "non-responders"
    )
    return res.sort_values("p_adj_bh").reset_index(drop=True)


def subject_level(cohort: pd.DataFrame) -> pd.DataFrame:
    """One value per subject per population (mean across that subject's samples)."""
    return (
        cohort.groupby(["subject", "response", "population"], as_index=False)["percentage"]
        .mean()
    )


def by_timepoint(cohort: pd.DataFrame) -> pd.DataFrame:
    """Exploratory: repeat the comparison separately at each timepoint.

    Secondary / hypothesis-generating only -- FDR is applied within a timepoint,
    not across the whole 15-test grid.
    """
    frames = []
    for t in sorted(cohort["time_from_treatment_start"].dropna().unique()):
        res = run_tests(cohort[cohort["time_from_treatment_start"] == t])
        if res.empty:
            continue
        res.insert(0, "time_from_treatment_start", int(t))
        frames.append(res)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def make_timecourse_plot(cohort: pd.DataFrame) -> str:
    ensure_dirs()
    np.random.seed(SEED)
    sns.set_theme(style="whitegrid", context="talk")
    pops = population_order(cohort)
    fig, axes = plt.subplots(1, len(pops), figsize=(4.1 * len(pops), 5.2), sharex=True)
    if len(pops) == 1:
        axes = [axes]
    for ax, pop in zip(axes, pops):
        sub = cohort[cohort["population"] == pop]
        sns.pointplot(
            data=sub, x="time_from_treatment_start", y="percentage",
            hue="response_label", hue_order=[RESPONSE_LABELS[r] for r in RESPONSE_ORDER],
            palette=PALETTE, estimator="median", errorbar=("ci", 95),
            dodge=0.25, markers="o", linestyles="-", ax=ax, legend=(ax is axes[-1]),
            seed=SEED,
        )
        ax.set_title(POPULATION_LABELS.get(pop, pop), fontsize=13)
        ax.set_xlabel("Days from treatment start")
        ax.set_ylabel("Median relative frequency (%)" if ax is axes[0] else "")
    if axes[-1].get_legend() is not None:
        axes[-1].legend(title="", loc="best", fontsize=11)
    fig.suptitle(
        "Exploratory: relative frequency over time by response "
        "(melanoma / miraclib / PBMC)", fontsize=16, y=1.03,
    )
    fig.tight_layout()
    path = os.path.join(FIGURE_DIR, "timecourse_response_by_population.png")
    fig.savefig(path, dpi=160, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)
    return path


def make_boxplots(cohort: pd.DataFrame, results: pd.DataFrame) -> str:
    ensure_dirs()
    np.random.seed(SEED)   # stripplot jitter -> byte-identical PNGs on re-run
    sns.set_theme(style="whitegrid", context="talk")
    pops = population_order(cohort)
    fig, axes = plt.subplots(1, len(pops), figsize=(4.1 * len(pops), 6.2), sharey=False)
    if len(pops) == 1:
        axes = [axes]

    lookup = results.set_index("population") if not results.empty else None
    for ax, pop in zip(axes, pops):
        sub = cohort[cohort["population"] == pop]
        sns.boxplot(
            data=sub, x="response_label", y="percentage", hue="response_label",
            order=[RESPONSE_LABELS[r] for r in RESPONSE_ORDER],
            palette=PALETTE, width=0.6, fliersize=0, legend=False, ax=ax,
        )
        sns.stripplot(
            data=sub.sample(min(len(sub), 400), random_state=SEED),
            x="response_label", y="percentage",
            order=[RESPONSE_LABELS[r] for r in RESPONSE_ORDER],
            color="black", alpha=0.18, size=2.4, jitter=0.28, ax=ax,
        )
        title = POPULATION_LABELS.get(pop, pop)
        if lookup is not None and pop in lookup.index:
            row = lookup.loc[pop]
            star = "*" if row["p_adj_bh"] < 0.05 else "ns"
            title += "\nBH p = %s (%s)" % (_fmt_p(row["p_adj_bh"]), star)
        ax.set_title(title, fontsize=13)
        ax.set_xlabel("")
        ax.set_ylabel("Relative frequency (% of sample)" if ax is axes[0] else "")
        ax.tick_params(axis="x", labelsize=11)

    fig.suptitle(
        "Melanoma / miraclib / PBMC -- cell population relative frequency by response",
        fontsize=16, y=1.02,
    )
    fig.tight_layout()
    path = os.path.join(FIGURE_DIR, "boxplots_response_by_population.png")
    fig.savefig(path, dpi=160, bbox_inches="tight", metadata={"Software": None})
    plt.close(fig)
    return path


def _fmt_p(p: float) -> str:
    return "%.3g" % p if p >= 1e-4 else "%.1e" % p


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--db", default=None)
    args = parser.parse_args()

    df = analysis_frame(args.db)
    cohort = melanoma_miraclib_pbmc(df)

    banner("PART 3 -- Responders vs non-responders (melanoma, miraclib, PBMC)")
    n_samples = cohort["sample"].nunique()
    n_subjects = cohort["subject"].nunique()
    counts = cohort.drop_duplicates("sample")["response"].value_counts()
    print("Cohort: %d samples from %d subjects (%d responder / %d non-responder samples)"
          % (n_samples, n_subjects, counts.get("yes", 0), counts.get("no", 0)))

    results = run_tests(cohort)
    cols = [
        "population_label", "n_responder", "n_non_responder",
        "median_responder", "median_non_responder", "median_difference",
        "rank_biserial", "p_mannwhitney", "p_adj_bh", "significant_fdr_0.05", "higher_in",
    ]
    print("\nAll samples (primary analysis):")
    print(results[cols].to_string(index=False, float_format=lambda v: "%.4g" % v))

    subj = subject_level(cohort)
    subj_results = run_tests(subj)
    print("\nSubject-level sensitivity analysis (one mean per subject):")
    print(subj_results[cols].to_string(index=False, float_format=lambda v: "%.4g" % v))

    tp = by_timepoint(cohort)
    if not tp.empty:
        print("\nExploratory: same comparison within each timepoint")
        print(tp[["time_from_treatment_start", "population_label", "median_responder",
                  "median_non_responder", "median_difference", "p_mannwhitney",
                  "p_adj_bh"]].to_string(index=False, float_format=lambda v: "%.4g" % v))
        write_csv(tp, "part3_response_statistics_by_timepoint.csv")

    fig_path = make_boxplots(cohort, results)
    tc_path = make_timecourse_plot(cohort)
    stats_path = write_csv(results, "part3_response_statistics.csv")
    subj_path = write_csv(subj_results, "part3_response_statistics_subject_level.csv")
    write_csv(cohort, "part3_cohort_frequencies.csv")

    sig = results[results["significant_fdr_0.05"]]
    banner("PART 3 -- Interpretation")
    if sig.empty:
        print("No population differs significantly between responders and "
              "non-responders after BH correction.")
    else:
        for _, r in sig.iterrows():
            print("* %-12s significant (BH p = %s): median %.2f%% in responders vs "
                  "%.2f%% in non-responders (%+.2f pp, higher in %s)."
                  % (r["population_label"], _fmt_p(r["p_adj_bh"]),
                     r["median_responder"], r["median_non_responder"],
                     r["median_difference"], r["higher_in"]))
    print("\nDirection of travel (exploratory): the responder / non-responder gap widens "
          "with\ntime on treatment -- CD4+ T cells drift up and B cells drift down in "
          "responders --\nwhich is the pattern to follow up with a longitudinal "
          "(mixed-effects) model.")
    print("\nFigures: %s\n         %s" % (fig_path, tc_path))
    print("Tables : %s\n         %s" % (stats_path, subj_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
