"""
avalanche_comparison.py
=======================
Comparative analysis of neuronal avalanches across brain states.

Supported comparisons
---------------------
1. sws  vs  non-sws  (other)
2. swsISR  vs  swsnonISR   (sws windows with / without InfraSlowRhythm)

For each condition pair the module produces:
  - P(S >= s) CCDF with power-law fit  (log-log)
  - Scalar metrics panel :
        tau_size,  tau_dur,
        gamma = (tau_dur - 1) / (tau_size - 1),
        tau_gamma_scaling  (expected from scaling relation),
        DCC = |tau_rapport - tau_gamma|
  - Statistical comparison :
        Z-test on tau_size / tau_dur between the two conditions
        KS two-sample test on raw size / duration distributions

Usage
-----
    from avalanche_comparison import AvalancheComparison

    # --- comparison 1 : sws vs other ---
    cmp = AvalancheComparison(session, bin_size=0.05, threshold=3)
    cmp.run_sws_vs_other(region="nr")
    fig = cmp.plot()
    fig.savefig("sws_vs_other.pdf")

    # --- comparison 2 : swsISR vs swsnonISR ---
    cmp2 = AvalancheComparison(session, bin_size=0.05, threshold=3)
    cmp2.run_swsISR_vs_swsnonISR(region="nr")
    fig2 = cmp2.plot()
    fig2.savefig("swsISR_vs_swsnonISR.pdf")
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.stats import ks_2samp
import powerlaw
from powerlaw import Fit
from scipy.optimize import OptimizeWarning

# ── try to import the project-level modules (adjust paths as needed) ──────────
try:
    import rg.data as rg_data          # noqa: F401  – used via rg.data.Regions
except ImportError:
    rg_data = None  # allow unit-testing without the full package


# ═══════════════════════════════════════════════════════════════════════════════
#  DATA CLASSES
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AvalancheResult:
    """Holds the full result for one brain-state condition."""
    label: str
    sizes:         np.ndarray            # raw avalanche sizes  (n,)
    durations:     np.ndarray            # raw avalanche durations (n,)
    # ── power-law fit on sizes ─────────────────────────────────────────────
    tau_size:      float = np.nan
    sigma_size:    float = np.nan
    xmin_size:     float = np.nan
    xmax_size:     float = np.nan
    R_size:        float = np.nan        # log-likelihood ratio vs exponential
    p_size:        float = np.nan
    ntail_size:    float = np.nan
    # ── power-law fit on durations ─────────────────────────────────────────
    tau_dur:       float = np.nan
    sigma_dur:     float = np.nan
    xmin_dur:      float = np.nan
    xmax_dur:      float = np.nan
    R_dur:         float = np.nan
    p_dur:         float = np.nan
    ntail_dur:     float = np.nan
    # ── derived criticality metrics ────────────────────────────────────────
    gamma_empirical:     float = np.nan   # (tau_dur-1)/(tau_size-1)
    tau_gamma_scaling:   float = np.nan   # from <S>~D^gamma  (direct regression)
    DCC:                 float = np.nan   # |tau_rapport - tau_gamma|

    # colour used in plots (set by AvalancheComparison)
    color: str = "C0"


@dataclass
class ComparisonStats:
    """Statistical comparison between two AvalancheResult objects."""
    # Z-tests (tau_A - tau_B) / sqrt(sigma_A^2 + sigma_B^2)
    z_tau_size:   float = np.nan
    p_z_tau_size: float = np.nan
    z_tau_dur:    float = np.nan
    p_z_tau_dur:  float = np.nan
    # KS two-sample
    ks_size:   float = np.nan
    p_ks_size: float = np.nan
    ks_dur:    float = np.nan
    p_ks_dur:  float = np.nan
    # difference in derived metrics
    delta_gamma:         float = np.nan
    delta_tau_scaling:   float = np.nan
    delta_DCC:           float = np.nan


# ═══════════════════════════════════════════════════════════════════════════════
#  LOW-LEVEL HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _powerlaw_fit(data: np.ndarray):
    """
    Fit a discrete power law via Clauset et al. 2009 (MLE).

    Returns
    -------
    tau, sigma, xmin, xmax, R_vs_exp, p_vs_exp, n_tail_proportion
    All NaN on failure.
    """
    data = np.asarray(data, dtype=float)
    data = data[data >= 1]

    nans = (np.nan,) * 7

    if len(np.unique(data)) < 2 or len(data) < 20:
        return nans

    data_int = data.astype(int)
    data_int = data_int[data_int >= 1]

    unique_vals = np.unique(data_int)
    if len(unique_vals) < 4:
        return nans

    most_common_freq = np.max(np.bincount(data_int)) / len(data_int)
    if most_common_freq > 0.90:
        return nans

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        warnings.simplefilter("ignore", UserWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            fit = Fit(data_int, discrete=True, verbose=False)
        except Exception:
            return nans

    tau   = fit.power_law.alpha
    sigma = fit.power_law.sigma
    xmin  = fit.xmin
    xmax  = fit.xmax if fit.xmax is not None else np.max(data_int)
    n_tail_prop = fit.n_tail / len(data)

    if np.isnan(tau):
        return nans

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        R, p = fit.distribution_compare(
            'power_law', 'exponential', normalized_ratio=True
        )

    return tau, sigma, xmin, xmax, R, p, n_tail_prop


def _ccdf(data: np.ndarray):
    """Compute empirical CCDF P(X >= x)."""
    data = np.sort(data)
    n    = len(data)
    ccdf = 1.0 - np.arange(1, n + 1) / n   # P(X >= x)  (right-continuous)
    return data, ccdf


def _z_test_tau(tau_a, sigma_a, tau_b, sigma_b):
    """
    Two-sided Z-test H0: tau_A == tau_B.
    sigma_a/b are standard errors from MLE.
    Returns (z, p_two_sided).
    """
    from scipy.stats import norm
    se = np.sqrt(sigma_a**2 + sigma_b**2)
    if se == 0 or np.isnan(se):
        return np.nan, np.nan
    z = (tau_a - tau_b) / se
    p = 2 * norm.sf(np.abs(z))
    return float(z), float(p)


def _gamma_from_scaling(sizes: np.ndarray, durations: np.ndarray) -> float:
    """
    Estimate gamma via OLS regression of log<S|D> on log D.
    Returns slope (= gamma), or NaN.
    """
    sizes     = np.asarray(sizes,     dtype=float)
    durations = np.asarray(durations, dtype=float)

    valid = (sizes >= 1) & (durations >= 1)
    S, D  = sizes[valid], durations[valid]

    unique_D = np.unique(D)
    if len(unique_D) < 4:
        return np.nan

    mean_S = np.array([S[D == d].mean() for d in unique_D])
    log_D  = np.log(unique_D)
    log_S  = np.log(mean_S)

    keep = np.isfinite(log_D) & np.isfinite(log_S)
    if keep.sum() < 4:
        return np.nan

    gamma, _ = np.polyfit(log_D[keep], log_S[keep], 1)
    return float(gamma)


def _fit_avalanche_result(
    label: str,
    sizes: np.ndarray,
    durations: np.ndarray,
    color: str = "C0",
) -> AvalancheResult:
    """Run all fits and derive criticality metrics for one condition."""
    res = AvalancheResult(
        label=label,
        sizes=np.asarray(sizes, dtype=float),
        durations=np.asarray(durations, dtype=float),
        color=color,
    )

    (res.tau_size, res.sigma_size, res.xmin_size,
     res.xmax_size, res.R_size, res.p_size, res.ntail_size) = _powerlaw_fit(sizes)

    (res.tau_dur, res.sigma_dur, res.xmin_dur,
     res.xmax_dur, res.R_dur, res.p_dur, res.ntail_dur) = _powerlaw_fit(durations)

    # gamma empirical from tau exponents
    if not (np.isnan(res.tau_size) or np.isnan(res.tau_dur)):
        denom = res.tau_size - 1
        if denom != 0:
            res.gamma_empirical = (res.tau_dur - 1) / denom

    # gamma from <S>~D^gamma scaling relation
    res.tau_gamma_scaling = _gamma_from_scaling(sizes, durations)

    # DCC = distance from criticality criterion
    if not (np.isnan(res.gamma_empirical) or np.isnan(res.tau_gamma_scaling)):
        res.DCC = abs(res.gamma_empirical - res.tau_gamma_scaling)

    return res


def _compare_two(a: AvalancheResult, b: AvalancheResult) -> ComparisonStats:
    """Compute all statistical comparisons between two conditions."""
    stats = ComparisonStats()

    # ── Z-tests on tau ────────────────────────────────────────────────────
    stats.z_tau_size, stats.p_z_tau_size = _z_test_tau(
        a.tau_size, a.sigma_size, b.tau_size, b.sigma_size
    )
    stats.z_tau_dur, stats.p_z_tau_dur = _z_test_tau(
        a.tau_dur, a.sigma_dur, b.tau_dur, b.sigma_dur
    )

    # ── KS two-sample tests ───────────────────────────────────────────────
    def _ks(x, y):
        x = np.asarray(x); y = np.asarray(y)
        x = x[np.isfinite(x) & (x >= 1)]
        y = y[np.isfinite(y) & (y >= 1)]
        if len(x) < 5 or len(y) < 5:
            return np.nan, np.nan
        res = ks_2samp(x, y)
        return float(res.statistic), float(res.pvalue)

    stats.ks_size,  stats.p_ks_size = _ks(a.sizes,     b.sizes)
    stats.ks_dur,   stats.p_ks_dur  = _ks(a.durations, b.durations)

    # ── differences in derived metrics ───────────────────────────────────
    stats.delta_gamma       = a.gamma_empirical   - b.gamma_empirical
    stats.delta_tau_scaling = a.tau_gamma_scaling - b.tau_gamma_scaling
    stats.delta_DCC         = a.DCC               - b.DCC

    return stats


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class AvalancheComparison:
    """
    High-level API for comparing avalanche statistics across brain states.

    Parameters
    ----------
    session       : neuroscience session object understood by rg.data.Regions
    bin_size      : float  — bin size (s) passed to get_avalanches
    threshold     : float  — threshold passed to get_avalanches
    get_avalanches: callable(FR, bin_size, threshold) -> (sizes, durations)
                    If None, the module will look for a globally available
                    `get_avalanches` function.
    """

    COLORS = {
        "sws":         "#2166ac",
        "other":       "#d6604d",
        "swsISR":      "#1a9850",
        "swsnonISR":   "#d73027",
    }

    def __init__(
        self,
        session,
        bin_size:       float = 0.05,
        threshold:      float = 3.0,
        get_avalanches=None,
    ):
        self.session      = session
        self.bin_size     = bin_size
        self.threshold    = threshold

        if get_avalanches is None:
            import builtins
            get_avalanches = getattr(builtins, "get_avalanches", None)
            if get_avalanches is None:
                # last resort: look in caller's global namespace
                import inspect
                frame = inspect.stack()[1][0]
                get_avalanches = frame.f_globals.get("get_avalanches")
        if get_avalanches is None:
            raise ValueError(
                "No get_avalanches function found. "
                "Pass it explicitly: AvalancheComparison(session, …, get_avalanches=fn)"
            )
        self._get_avalanches = get_avalanches

        self.results:  list[AvalancheResult] = []
        self.stats:    Optional[ComparisonStats] = None
        self._title:   str = ""

    # ── spike / FR helpers ────────────────────────────────────────────────────

    def _get_FR(self, region: str, state: str, shift: bool = True,
                events: Optional[str] = None):
        """
        Return firing-rate matrix for `region` restricted to `state`.
        Uses rg.data.Regions (with optional ISR events).
        """
        import rg.data as rg_data  # local import to keep the module importable

        states_needed = ["sws", "rem"] if state in ("sws", "rem") else [state, "sws"]
        # deduplicate while preserving order
        seen = set(); states_needed_u = []
        for s in states_needed:
            if s not in seen:
                seen.add(s); states_needed_u.append(s)

        kwargs = dict(states=states_needed_u)
        if events is not None:
            kwargs["events"] = events

        R = rg_data.Regions(self.session, **kwargs)
        FR = R.spikes(regs=[region], state=state, shift=shift)
        return FR

    def _avalanches_for(self, FR):
        sizes, durations = self._get_avalanches(FR, self.bin_size, self.threshold)
        return np.asarray(sizes, dtype=float), np.asarray(durations, dtype=float)

    # ── public run methods ────────────────────────────────────────────────────

    def run_sws_vs_other(self, region: str = "nr"):
        """
        Compare SWS avalanches vs non-SWS (other) avalanches.
        """
        self._title = f"SWS vs Other — [{region}]"

        # SWS
        FR_sws = self._get_FR(region, state="sws")
        s_sws, d_sws = self._avalanches_for(FR_sws)

        # Other (non-sws) — use state='other' if supported, else 'wake'
        try:
            FR_other = self._get_FR(region, state="other")
        except Exception:
            FR_other = self._get_FR(region, state="wake")
        s_other, d_other = self._avalanches_for(FR_other)

        self.results = [
            _fit_avalanche_result("SWS",   s_sws,   d_sws,   self.COLORS["sws"]),
            _fit_avalanche_result("Other", s_other, d_other, self.COLORS["other"]),
        ]
        self.stats = _compare_two(*self.results)
        return self

    def run_swsISR_vs_swsnonISR(self, region: str = "nr"):
        """
        Compare SWS-ISR avalanches vs SWS-non-ISR avalanches.
        nonISR = sws windows EXCLUDING InfraSlowRhythm events.
        """
        self._title = f"SWS-ISR vs SWS-nonISR — [{region}]"

        # ISR windows: Regions with events="InfraSlowRythm"  (typo in original API)
        FR_isr = self._get_FR(region, state="sws", events="InfraSlowRythm")
        s_isr, d_isr = self._avalanches_for(FR_isr)

        # non-ISR: sws WITHOUT events keyword (full sws), then subtract ISR
        # Some rg implementations accept a negation flag; we handle both cases.
        try:
            FR_nonisr = self._get_FR(region, state="sws", events="~InfraSlowRythm")
            s_nonisr, d_nonisr = self._avalanches_for(FR_nonisr)
        except Exception:
            # Fallback: get full sws and remove ISR-tagged avalanches is not
            # trivially possible at the spike level → just use full sws as proxy
            # and warn the user.
            warnings.warn(
                "Could not obtain sws\\ISR via events='~InfraSlowRythm'. "
                "Using full-sws as swsnonISR proxy. "
                "Override _get_FR_swsnonISR() for a clean exclusion.",
                UserWarning, stacklevel=2,
            )
            FR_nonisr = self._get_FR(region, state="sws")
            s_nonisr, d_nonisr = self._avalanches_for(FR_nonisr)

        self.results = [
            _fit_avalanche_result("SWS-ISR",    s_isr,    d_isr,    self.COLORS["swsISR"]),
            _fit_avalanche_result("SWS-nonISR", s_nonisr, d_nonisr, self.COLORS["swsnonISR"]),
        ]
        self.stats = _compare_two(*self.results)
        return self

    # ── plotting ──────────────────────────────────────────────────────────────

    def plot(self, figsize=(16, 10)) -> plt.Figure:
        """
        Generate the full comparison figure.

        Layout
        ------
        Row 0 :  CCDF size (A)  |  CCDF duration (B)
        Row 1 :  Metrics bar chart (C)   |   Stats table (D)
        """
        if not self.results:
            raise RuntimeError("Call run_*() before plot().")

        fig = plt.figure(figsize=figsize, constrained_layout=True)
        fig.suptitle(self._title, fontsize=14, fontweight="bold", y=1.01)

        gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)
        ax_size  = fig.add_subplot(gs[0, 0])
        ax_dur   = fig.add_subplot(gs[0, 1])
        ax_met   = fig.add_subplot(gs[1, 0])
        ax_stat  = fig.add_subplot(gs[1, 1])

        self._plot_ccdf(ax_size, which="size")
        self._plot_ccdf(ax_dur,  which="duration")
        self._plot_metrics(ax_met)
        self._plot_stats_table(ax_stat)

        return fig

    # ── internal plot helpers ─────────────────────────────────────────────────

    def _plot_ccdf(self, ax: plt.Axes, which: str = "size"):
        """Log-log CCDF + power-law fit for sizes or durations."""
        is_size = (which == "size")
        for res in self.results:
            raw  = res.sizes     if is_size else res.durations
            tau  = res.tau_size  if is_size else res.tau_dur
            xmin = res.xmin_size if is_size else res.xmin_dur
            xmax = res.xmax_size if is_size else res.xmax_dur

            raw = raw[np.isfinite(raw) & (raw >= 1)]
            if len(raw) < 5:
                continue

            x_ccdf, y_ccdf = _ccdf(raw)
            ax.plot(x_ccdf, y_ccdf, ".", color=res.color,
                    alpha=0.4, markersize=3, rasterized=True)

            # ── power-law line in the fitted range ──────────────────────
            if np.isfinite(tau) and np.isfinite(xmin):
                xs  = np.logspace(np.log10(xmin), np.log10(xmax), 200)
                # normalise so that the line passes through the CCDF at xmin
                idx   = np.searchsorted(x_ccdf, xmin)
                y0    = y_ccdf[min(idx, len(y_ccdf) - 1)]
                slope = -(tau - 1)
                ys    = y0 * (xs / xmin) ** slope
                label = f"{res.label}  τ={tau:.2f}±{(res.sigma_size if is_size else res.sigma_dur):.2f}"
                ax.plot(xs, ys, "-", color=res.color, linewidth=2, label=label)

        xlabel = "Avalanche size $s$" if is_size else "Avalanche duration $d$"
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(r"$P(S \geq s)$" if is_size else r"$P(D \geq d)$", fontsize=11)
        ax.set_title("Size CCDF" if is_size else "Duration CCDF", fontsize=12)
        ax.legend(fontsize=9, framealpha=0.7)
        ax.grid(True, which="both", alpha=0.25, linestyle="--")

    def _plot_metrics(self, ax: plt.Axes):
        """Grouped bar chart of criticality metrics per condition."""
        metrics = [
            ("γ emp.\n(τ_dur-1)/(τ_size-1)", "gamma_empirical"),
            ("γ scaling\n⟨S⟩~D^γ",            "tau_gamma_scaling"),
            ("DCC\n|γ_emp - γ_scal|",          "DCC"),
        ]
        n_metrics = len(metrics)
        n_cond    = len(self.results)
        width     = 0.8 / n_cond
        x         = np.arange(n_metrics)

        for i, res in enumerate(self.results):
            vals = [getattr(res, attr) for _, attr in metrics]
            offset = (i - (n_cond - 1) / 2) * width
            bars = ax.bar(x + offset, vals, width * 0.9,
                          color=res.color, label=res.label, alpha=0.85,
                          edgecolor="white", linewidth=0.5)
            for bar, v in zip(bars, vals):
                if np.isfinite(v):
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_height() + 0.01,
                            f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)

        ax.set_xticks(x)
        ax.set_xticklabels([m for m, _ in metrics], fontsize=9.5)
        ax.set_ylabel("Value", fontsize=11)
        ax.set_title("Criticality metrics", fontsize=12)
        ax.legend(fontsize=9, framealpha=0.7)
        ax.axhline(0, color="black", linewidth=0.6)
        ax.grid(axis="y", alpha=0.25, linestyle="--")

    def _plot_stats_table(self, ax: plt.Axes):
        """Table summarising all statistical tests."""
        ax.axis("off")
        if self.stats is None:
            return

        s  = self.stats
        la = self.results[0].label
        lb = self.results[1].label

        def _fmt(v):
            return f"{v:.4f}" if np.isfinite(v) else "—"

        def _star(p):
            if not np.isfinite(p): return ""
            if p < 0.001: return "***"
            if p < 0.01:  return "**"
            if p < 0.05:  return "*"
            return "n.s."

        rows = [
            ["Test", "Statistic", "p-value", "Sig."],
            # ── Z-tests ──────────────────────────────────────────────
            [f"Z-test τ_size\n({la} vs {lb})",
             f"z={_fmt(s.z_tau_size)}", _fmt(s.p_z_tau_size), _star(s.p_z_tau_size)],
            [f"Z-test τ_dur\n({la} vs {lb})",
             f"z={_fmt(s.z_tau_dur)}", _fmt(s.p_z_tau_dur), _star(s.p_z_tau_dur)],
            # ── KS tests ─────────────────────────────────────────────
            [f"KS sizes\n({la} vs {lb})",
             f"D={_fmt(s.ks_size)}", _fmt(s.p_ks_size), _star(s.p_ks_size)],
            [f"KS durations\n({la} vs {lb})",
             f"D={_fmt(s.ks_dur)}", _fmt(s.p_ks_dur), _star(s.p_ks_dur)],
            # ── metric diffs ─────────────────────────────────────────
            [f"Δγ empirical\n({la}−{lb})",
             f"{_fmt(s.delta_gamma)}", "—", ""],
            [f"Δγ scaling\n({la}−{lb})",
             f"{_fmt(s.delta_tau_scaling)}", "—", ""],
            [f"ΔDCC\n({la}−{lb})",
             f"{_fmt(s.delta_DCC)}", "—", ""],
        ]

        col_widths = [0.38, 0.24, 0.20, 0.10]
        n_rows     = len(rows)
        row_h      = 1.0 / n_rows

        for r_idx, row in enumerate(rows):
            y = 1.0 - (r_idx + 0.5) * row_h
            is_header = (r_idx == 0)
            x_cursor  = 0.0

            # row background
            bg = "#e8e8e8" if is_header else ("#f7f7f7" if r_idx % 2 == 0 else "white")
            rect = plt.Rectangle((0, 1.0 - (r_idx + 1) * row_h),
                                  1.0, row_h,
                                  transform=ax.transAxes,
                                  color=bg, zorder=0, clip_on=False)
            ax.add_patch(rect)

            for c_idx, (cell, cw) in enumerate(zip(row, col_widths)):
                xc = x_cursor + cw / 2
                fw = "bold" if is_header else "normal"
                fs = 8.5 if is_header else 8

                # colour p-value column based on significance
                fc = "black"
                if c_idx == 2 and not is_header:
                    try:
                        pv = float(cell)
                        if pv < 0.001:   fc = "#006d2c"
                        elif pv < 0.01:  fc = "#31a354"
                        elif pv < 0.05:  fc = "#74c476"
                        else:            fc = "#d73027"
                    except ValueError:
                        pass

                ax.text(xc, y, cell, ha="center", va="center",
                        fontsize=fs, fontweight=fw, color=fc,
                        transform=ax.transAxes, wrap=True,
                        multialignment="center")
                x_cursor += cw

        # outer border
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title("Statistical comparison", fontsize=12, pad=6)

    # ── convenience: print summary ────────────────────────────────────────────

    def summary(self):
        """Print a text summary of all metrics and comparisons."""
        print(f"\n{'='*60}")
        print(f"  {self._title}")
        print(f"{'='*60}")
        for res in self.results:
            print(f"\n  ── {res.label} (n={len(res.sizes)}) ──")
            print(f"     τ_size          = {res.tau_size:.4f} ± {res.sigma_size:.4f}"
                  f"   (xmin={res.xmin_size}, p_vs_exp={res.p_size:.3f})")
            print(f"     τ_dur           = {res.tau_dur:.4f} ± {res.sigma_dur:.4f}"
                  f"   (xmin={res.xmin_dur}, p_vs_exp={res.p_dur:.3f})")
            print(f"     γ empirical     = {res.gamma_empirical:.4f}")
            print(f"     γ scaling       = {res.tau_gamma_scaling:.4f}")
            print(f"     DCC             = {res.DCC:.4f}")

        if self.stats:
            s = self.stats
            la, lb = self.results[0].label, self.results[1].label
            print(f"\n  ── Statistical tests ({la} vs {lb}) ──")
            print(f"     Z τ_size  z={s.z_tau_size:.3f}  p={s.p_z_tau_size:.4f}")
            print(f"     Z τ_dur   z={s.z_tau_dur:.3f}  p={s.p_z_tau_dur:.4f}")
            print(f"     KS sizes  D={s.ks_size:.3f}   p={s.p_ks_size:.4f}")
            print(f"     KS durs   D={s.ks_dur:.3f}    p={s.p_ks_dur:.4f}")
        print()


# ═══════════════════════════════════════════════════════════════════════════════
#  BATCH HELPER  – run both comparisons for a list of sessions
# ═══════════════════════════════════════════════════════════════════════════════

def run_all_comparisons(
    sessions: list,
    bin_size:       float = 0.05,
    threshold:      float = 3.0,
    region:         str   = "nr",
    get_avalanches=None,
    save_dir:       Optional[str] = None,
) -> dict[str, list]:
    """
    Run both sws-vs-other and swsISR-vs-swsnonISR for every session.

    Returns
    -------
    dict with keys "sws_vs_other" and "swsISR_vs_swsnonISR",
    each mapping to a list of AvalancheComparison objects.

    If save_dir is given, figures are saved as PDF there.
    """
    import os

    out = {"sws_vs_other": [], "swsISR_vs_swsnonISR": []}

    for i, sess in enumerate(sessions):
        print(f"[{i+1}/{len(sessions)}] {sess}")

        for key, run_fn in [
            ("sws_vs_other",       lambda c: c.run_sws_vs_other(region)),
            ("swsISR_vs_swsnonISR", lambda c: c.run_swsISR_vs_swsnonISR(region)),
        ]:
            try:
                cmp = AvalancheComparison(
                    sess, bin_size=bin_size, threshold=threshold,
                    get_avalanches=get_avalanches,
                )
                run_fn(cmp)
                out[key].append(cmp)

                if save_dir is not None:
                    os.makedirs(save_dir, exist_ok=True)
                    fname = os.path.join(save_dir, f"sess{i:03d}_{key}.pdf")
                    fig   = cmp.plot()
                    fig.savefig(fname, bbox_inches="tight")
                    plt.close(fig)
                    print(f"   saved → {fname}")

            except Exception as exc:
                warnings.warn(f"Session {i} [{key}] failed: {exc}", RuntimeWarning)

    return out



# SWS vs Other
cmp = AvalancheComparison(session, bin_size=0.05, threshold=3.0,
                           get_avalanches=get_avalanches)
cmp.run_sws_vs_other(region="nr")
fig = cmp.plot()

# SWS-ISR vs SWS-nonISR
cmp2 = AvalancheComparison(session, bin_size=0.05, threshold=3.0,
                            get_avalanches=get_avalanches)
cmp2.run_swsISR_vs_swsnonISR(region="nr")
fig2 = cmp2.plot()