import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Ellipse
from adjustText import adjust_text

# =========================
# helper: simple gaussian KDE
# =========================
def get_kde_curve(x, gridsize=300, extend=0.08):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) == 0:
        return np.array([0, 1]), np.array([0, 0])

    if len(x) == 1:
        return np.array([x[0] - 1e-3, x[0] + 1e-3]), np.array([1, 1])

    xmin, xmax = x.min(), x.max()
    xrng = xmax - xmin
    if xrng == 0:
        xrng = max(abs(xmin), 1.0) * 0.1

    grid = np.linspace(xmin - extend * xrng, xmax + extend * xrng, gridsize)

    # Silverman's rule
    std = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    sigma = min(std, iqr / 1.34) if iqr > 0 else std
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = max(std, xrng / 20, 1e-3)

    bw = 0.9 * sigma * (len(x) ** (-1 / 5))
    if not np.isfinite(bw) or bw <= 0:
        bw = max(xrng / 50, 1e-3)

    diff = (grid[:, None] - x[None, :]) / bw
    density = np.exp(-0.5 * diff**2).sum(axis=1) / (len(x) * bw * np.sqrt(2 * np.pi))

    return grid, density


# =========================
# helper: draw one ellipse around one group
# =========================
def add_group_ellipse(ax, df, xcol, ycol, edgecolor, pad_frac=0.12, min_size_frac=0.06, lw=2.2):
    """
    Draw one large ellipse around all points in df.
    pad_frac: padding relative to group range
    min_size_frac: minimum ellipse size relative to full axis range
    """
    if df.empty:
        return

    x = df[xcol].values
    y = df[ycol].values

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    xspan = xmax - xmin
    yspan = ymax - ymin

    # Use full data range to ensure tiny groups still get visible ellipse
    ax_xlim = ax.get_xlim()
    ax_ylim = ax.get_ylim()
    full_xspan = ax_xlim[1] - ax_xlim[0]
    full_yspan = ax_ylim[1] - ax_ylim[0]

    min_w = full_xspan * min_size_frac
    min_h = full_yspan * min_size_frac

    pad_x = max(xspan * pad_frac, min_w / 2)
    pad_y = max(yspan * pad_frac, min_h / 2)

    width = max(xspan + 2 * pad_x, min_w)
    height = max(yspan + 2 * pad_y, min_h)

    center_x = (xmin + xmax) / 2
    center_y = (ymin + ymax) / 2

    ell = Ellipse(
        (center_x, center_y),
        width=width,
        height=height,
        angle=0,
        facecolor="none",
        edgecolor=edgecolor,
        linewidth=lw,
        linestyle="-",
        alpha=0.95,
        zorder=2.8
    )
    ax.add_patch(ell)


# =========================
# paths
# =========================
task_id = int(os.getenv("SLURM_ARRAY_TASK_ID", "1"))

if task_id == 1:
    plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula/NMF"
    table_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF"
    data_suffix = ""
else:
    plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula_extracellular/NMF"
    table_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF_extracellular"
    data_suffix = ""

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(table_dir, exist_ok=True)

# =========================
# read saved loadings
# =========================
loadings_file = os.path.join(table_dir, f"NMF_H_loadings{data_suffix}.csv")
lr_loadings = pd.read_csv(loadings_file, index_col=0)

print("Loaded:", loadings_file)
print("Available columns:", lr_loadings.columns.tolist())

# =========================
# settings
# =========================
factor_x = "Factor3"
factor_y = "Factor4"
top_n = 20

if factor_x not in lr_loadings.columns or factor_y not in lr_loadings.columns:
    raise ValueError(
        f"{factor_x} or {factor_y} not found. Available columns: {lr_loadings.columns.tolist()}"
    )

# =========================
# prepare dataframe
# =========================
scatter_df = (
    lr_loadings[[factor_x, factor_y]]
    .copy()
    .reset_index()
    .rename(columns={"index": "lr_pair"})
)

top_x = set(lr_loadings[factor_x].nlargest(top_n).index)
top_y = set(lr_loadings[factor_y].nlargest(top_n).index)
highlight_pairs = top_x.union(top_y)

def classify_pair(x):
    in_x = x in top_x
    in_y = x in top_y
    if in_x and in_y:
        return "Top in both"
    elif in_x:
        return f"Top in {factor_x}"
    elif in_y:
        return f"Top in {factor_y}"
    else:
        return "Other"

scatter_df["group"] = scatter_df["lr_pair"].apply(classify_pair)
scatter_df["label_score"] = scatter_df[factor_x] + scatter_df[factor_y]

# save highlighted table
highlight_df = scatter_df[scatter_df["lr_pair"].isin(highlight_pairs)].copy()
highlight_df[f"in_{factor_x}_top{top_n}"] = highlight_df["lr_pair"].isin(top_x)
highlight_df[f"in_{factor_y}_top{top_n}"] = highlight_df["lr_pair"].isin(top_y)
highlight_df = highlight_df.sort_values("label_score", ascending=False)

highlight_table = os.path.join(
    table_dir, f"{factor_x}_{factor_y}_top{top_n}_highlighted_pairs.csv"
)
highlight_df.to_csv(highlight_table, index=False)
print("Saved:", highlight_table)

# =========================
# plot style
# =========================
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.linewidth": 1.0
})

# =========================
# layout
# =========================
fig = plt.figure(figsize=(13.5, 10.5))
gs = fig.add_gridspec(
    2, 2,
    width_ratios=(5.4, 0.8),
    height_ratios=(0.8, 5.4),
    wspace=0.05,
    hspace=0.05
)

ax_top = fig.add_subplot(gs[0, 0])
ax = fig.add_subplot(gs[1, 0], sharex=ax_top)
ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

fig.patch.set_facecolor("white")
ax.set_facecolor("white")
ax_top.set_facecolor("white")
ax_right.set_facecolor("white")

fig.suptitle(
    f"{factor_x} vs {factor_y}: Top {top_n} ligand–receptor pairs",
    fontsize=22,
    y=0.98
)

# light grid only on main scatter
ax.grid(True, which="major", color="#EAEAEA", linewidth=1.0)
ax.set_axisbelow(True)

# =========================
# colors and sizes
# =========================
color_map = {
    f"Top in {factor_x}": "#2C7FB8",
    f"Top in {factor_y}": "#D95F5F",
    "Top in both": "#2CA25F"
}

size_map = {
    "Other": 26,
    f"Top in {factor_x}": 90,
    f"Top in {factor_y}": 90,
    "Top in both": 120
}

highlight_order = [f"Top in {factor_x}", f"Top in {factor_y}", "Top in both"]

# =========================
# main scatter
# =========================
base_df = scatter_df[scatter_df["group"] == "Other"]

ax.scatter(
    base_df[factor_x],
    base_df[factor_y],
    s=size_map["Other"],
    color="#CFCFCF",
    alpha=0.40,
    edgecolors="none",
    zorder=1
)

for grp in highlight_order:
    sub = scatter_df[scatter_df["group"] == grp]
    if sub.empty:
        continue

    ax.scatter(
        sub[factor_x],
        sub[factor_y],
        s=size_map[grp],
        color=color_map[grp],
        alpha=0.95,
        edgecolors="white",
        linewidths=0.8,
        label=grp,
        zorder=4
    )

# reference lines
ax.axhline(0, color="#BDBDBD", linewidth=1.0, linestyle="--", zorder=0)
ax.axvline(0, color="#BDBDBD", linewidth=1.0, linestyle="--", zorder=0)

# set limits first so ellipse helper can use them
xvals = scatter_df[factor_x].values
yvals = scatter_df[factor_y].values
xpad = 0.06 * (xvals.max() - xvals.min() if xvals.max() > xvals.min() else 1)
ypad = 0.06 * (yvals.max() - yvals.min() if yvals.max() > yvals.min() else 1)
ax.set_xlim(xvals.min() - xpad, xvals.max() + xpad)
ax.set_ylim(yvals.min() - ypad, yvals.max() + ypad)

# =========================
# draw ONLY 3 big ellipses
# =========================
for grp in highlight_order:
    sub = scatter_df[scatter_df["group"] == grp]
    add_group_ellipse(
        ax=ax,
        df=sub,
        xcol=factor_x,
        ycol=factor_y,
        edgecolor=color_map[grp],
        pad_frac=0.18,
        min_size_frac=0.08,
        lw=2.2
    )

# =========================
# overall smooth density curves only
# =========================
grid_x, dens_x = get_kde_curve(scatter_df[factor_x].values)
ax_top.fill_between(grid_x, dens_x, color="#D9D9D9", alpha=0.75, linewidth=0)
ax_top.plot(grid_x, dens_x, color="#AFAFAF", linewidth=1.6)

grid_y, dens_y = get_kde_curve(scatter_df[factor_y].values)
ax_right.fill_betweenx(grid_y, 0, dens_y, color="#D9D9D9", alpha=0.75, linewidth=0)
ax_right.plot(dens_y, grid_y, color="#AFAFAF", linewidth=1.6)

# =========================
# labels
# =========================
label_df = scatter_df[scatter_df["lr_pair"].isin(highlight_pairs)].copy()
label_df = label_df.sort_values("label_score", ascending=False)

texts = []
for _, row in label_df.iterrows():
    txt = ax.text(
        row[factor_x],
        row[factor_y],
        row["lr_pair"],
        fontsize=9,
        color="#222222",
        zorder=6
    )
    txt.set_path_effects([
        pe.withStroke(linewidth=2.2, foreground="white", alpha=0.9)
    ])
    texts.append(txt)

adjust_text(
    texts,
    ax=ax,
    expand=(1.15, 1.25),
    force_text=(0.8, 1.0),
    force_static=(0.4, 0.6),
    force_pull=(0.1, 0.1),
    arrowprops=dict(
        arrowstyle="-",
        color="#888888",
        lw=0.6,
        alpha=0.7
    )
)

# =========================
# axes and legend
# =========================
ax.set_xlabel(f"{factor_x} loading")
ax.set_ylabel(f"{factor_y} loading")
ax.set_title(f"{factor_x} vs {factor_y}: Top {top_n} ligand–receptor pairs", pad=12)

# main axis
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#888888")
ax.spines["bottom"].set_color("#888888")

# top density axis
ax_top.spines["top"].set_visible(False)
ax_top.spines["right"].set_visible(False)
ax_top.spines["left"].set_visible(False)
ax_top.spines["bottom"].set_color("#BBBBBB")
ax_top.tick_params(axis="x", labelbottom=False, bottom=False)
ax_top.tick_params(axis="y", left=False, labelleft=False)

# right density axis
ax_right.spines["top"].set_visible(False)
ax_right.spines["right"].set_visible(False)
ax_right.spines["bottom"].set_visible(False)
ax_right.spines["left"].set_color("#BBBBBB")
ax_right.tick_params(axis="y", labelleft=False, left=False)
ax_right.tick_params(axis="x", bottom=False, labelbottom=False)

legend = ax.legend(
    title=None,
    frameon=True,
    facecolor="white",
    edgecolor="none",
    framealpha=0.85,
    loc="upper right",
    fontsize=11,
    markerscale=0.9,
    handletextpad=0.4,
    borderpad=0.3,
    labelspacing=0.4
)

plt.tight_layout(rect=[0, 0, 1, 0.95])

pdf_file = os.path.join(
    plot_dir, f"{factor_x}_{factor_y}_top{top_n}_scatter_density_3ellipses{data_suffix}.pdf"
)
png_file = os.path.join(
    plot_dir, f"{factor_x}_{factor_y}_top{top_n}_scatter_density_3ellipses{data_suffix}.png"
)

plt.savefig(pdf_file, dpi=300, bbox_inches="tight")
plt.savefig(png_file, dpi=300, bbox_inches="tight")
plt.close()

print("✅ Saved:", pdf_file)
print("✅ Saved:", png_file)