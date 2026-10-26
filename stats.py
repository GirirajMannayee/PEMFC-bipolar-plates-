"""Statistics reproduced from summary data (mean, SD, n) as in manuscript Sec. 2.9.

All routines work from tabulated mean +/- SD because raw replicates are not
distributed with the manuscript. With raw replicates, use the *_raw variants
(statsmodels / scipy) - the summary-statistics results are mathematically
identical for balanced designs.
"""
from __future__ import annotations
import itertools
import numpy as np
import pandas as pd
from scipy import stats


def oneway_from_summary(means, sds, n):
    """One-way ANOVA from group means/SDs (equal n). Returns dict."""
    m, s = np.asarray(means, float), np.asarray(sds, float)
    k = len(m)
    grand = m.mean()
    ssb = n * np.sum((m - grand) ** 2)
    ssw = (n - 1) * np.sum(s ** 2)
    dfb, dfw = k - 1, k * (n - 1)
    msb, msw = ssb / dfb, ssw / dfw
    F = msb / msw
    return dict(F=F, p=stats.f.sf(F, dfb, dfw), eta2=ssb / (ssb + ssw),
                dfb=dfb, dfw=dfw, mse=msw)


def tukey_from_summary(names, means, sds, n, alpha=0.05):
    """Tukey HSD from summary data. Returns (pairwise DataFrame, HSD, letters dict)."""
    m = np.asarray(means, float)
    res = oneway_from_summary(means, sds, n)
    k, dfw, mse = len(m), res["dfw"], res["mse"]
    q = stats.studentized_range.ppf(1 - alpha, k, dfw)
    hsd = q * np.sqrt(mse / n)
    rows = []
    for (i, a), (j, b) in itertools.combinations(list(enumerate(names)), 2):
        diff = m[i] - m[j]
        qstat = abs(diff) / np.sqrt(mse / n)
        rows.append(dict(a=a, b=b, diff=diff, q=qstat,
                         p=stats.studentized_range.sf(qstat, k, dfw),
                         significant=abs(diff) > hsd))
    return pd.DataFrame(rows), hsd, compact_letters(names, m, hsd)


def compact_letters(names, means, hsd):
    """Compact letter display; 'a' = highest mean. Simple insert/absorb algorithm."""
    order = np.argsort(-np.asarray(means))
    nm = [names[i] for i in order]
    mv = np.asarray(means)[order]
    groups = []  # list of sets of names that are mutually non-significant
    for i in range(len(nm)):
        grp = {nm[i]}
        for j in range(i + 1, len(nm)):
            if mv[i] - mv[j] <= hsd:
                grp.add(nm[j])
        groups.append(grp)
    # remove groups that are subsets of another group
    uniq = []
    for g in groups:
        if not any(g < h for h in groups) and g not in uniq:
            uniq.append(g)
    letters = {n: "" for n in names}
    for li, g in enumerate(uniq):
        for n in g:
            letters[n] += "abcdefgh"[li]
    return letters


def holm(pvals):
    """Holm step-down adjusted p-values (same order as input)."""
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        adj[idx] = min(1.0, running)
    return adj


def factorial_2x2(means, sds, n, log=False):
    """2x2 factorial (Cu x CNT) contrasts from summary data.

    means/sds ordered [S1 (none), S2 (Cu), S3 (CNT), S4 (Cu+CNT)].
    log=True: analyse ln(response); SD(ln) ~= SD/mean (delta method), effects
    returned as ratios. Contrast SEs use the pooled error variance (4*(n-1) df).
    """
    m, s = np.asarray(means, float), np.asarray(sds, float)
    if log:
        s = s / m
        m = np.log(m)
    df = 4 * (n - 1)
    mse = np.mean(s ** 2)
    tcrit = stats.t.ppf(0.975, df)
    contrasts = {
        "cu": (np.array([-1, 1, -1, 1]) / 2, np.sqrt(mse / n)),
        "cnt": (np.array([-1, -1, 1, 1]) / 2, np.sqrt(mse / n)),
        "interaction": (np.array([1, -1, -1, 1]), np.sqrt(4 * mse / n)),
    }
    out = {}
    for name, (c, se) in contrasts.items():
        est = float(c @ m)
        lo, hi = est - tcrit * se, est + tcrit * se
        t = est / se
        p = 2 * stats.t.sf(abs(t), df)
        if log:
            est, lo, hi = np.exp([est, lo, hi])
        out[name] = dict(estimate=est, lo=lo, hi=hi, t=t, p=p)
    ssb = n * np.sum((m - m.mean()) ** 2)
    out["eta2"] = ssb / (ssb + (n - 1) * np.sum(s ** 2))
    out["df"] = df
    return out
