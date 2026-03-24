import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from adjustText import adjust_text

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

fig, ax = plt.subplots(figsize=(9, 7.5))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# light grid
ax.grid(True, which="major", color="#EAEAEA", linewidth=1.0)
ax.set_axisbelow(True)

# =========================
# background points
# =========================
base_df = scatter_df[scatter_df["group"] == "Other"]
ax.scatter(
    base_df[factor_x],
    base_df[factor_y],
    s=26,
    color="#CFCFCF",
    alpha=0.45,
    edgecolors="none",
    zorder=1
)

# =========================
# highlighted points
# =========================
color_map = {
    f"Top in {factor_x}": "#2C7FB8",
    f"Top in {factor_y}": "#D95F5F",
    "Top in both": "#2CA25F"
}

size_map = {
    f"Top in {factor_x}": 90,
    f"Top in {factor_y}": 90,
    "Top in both": 120
}

for grp in [f"Top in {factor_x}", f"Top in {factor_y}", "Top in both"]:
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
        zorder=3
    )

# subtle reference lines
ax.axhline(0, color="#BDBDBD", linewidth=1.0, linestyle="--", zorder=0)
ax.axvline(0, color="#BDBDBD", linewidth=1.0, linestyle="--", zorder=0)

# =========================
# labels (no white bbox)
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
        zorder=4
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

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#888888")
ax.spines["bottom"].set_color("#888888")

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

plt.tight_layout()

pdf_file = os.path.join(
    plot_dir, f"{factor_x}_{factor_y}_top{top_n}_scatter_pretty{data_suffix}.pdf"
)
png_file = os.path.join(
    plot_dir, f"{factor_x}_{factor_y}_top{top_n}_scatter_pretty{data_suffix}.png"
)

plt.savefig(pdf_file, dpi=300, bbox_inches="tight")
plt.savefig(png_file, dpi=300, bbox_inches="tight")
plt.close()

print("✅ Saved:", pdf_file)
print("✅ Saved:", png_file)