# -*- coding: utf-8 -*-
"""
Generate a publication-ready non-normalized top-12 stacked abundance bar plot.

Dataset
-------
KEGG level-3 functional abundance matrix

Run in Spyder
-------------
1. Install dependencies listed in repository-root requirements.txt.
2. Run from any working directory (relative input paths resolve beside script).
3. Outputs will be written to repository-root output/figures/.

Main plotting rule
------------------
This script selects the top 12 displayed features by mean relative abundance
after label cleaning. Existing source rows such as "Others" are retained if
they are among the selected top 12. Lower-abundance features outside the
displayed set are NOT pooled back into "Others", and the plotted values are
NOT normalized to 100%.

Interpretation
--------------
The bars show the abundance contributed by the displayed top-12 features within
the total annotated profile. Therefore, bars may end below 100%. This is
intentional and avoids inflating "Others" or the displayed features.
"""

from pathlib import Path
import re
import textwrap
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager


# =============================================================================
# User settings
# =============================================================================

INPUT_CSV = r"../data/functional/kegg_l3.csv"
OUTPUT_DIR = "../output/figures"

FIGURE_PREFIX = "Figure_2_KEGG_Level3"
FEATURE_COLUMN = "KEGG_Level3"
LEGEND_TITLE = "KEGG level 3 pathway"
PLOT_TITLE = "KEGG level 3 functional abundance"

TOP_N = 12
SHOW_TITLE = False

# Keep this False for the main manuscript figure unless you specifically want
# all labels containing "Unclassified" to collapse into one "Unresolved" row.
POOL_UNCLASSIFIED_AS_UNRESOLVED = False

# If True, only exact terms such as "Unresolved", "Unassigned", "Unannotated",
# and "Unknown" are standardised to "Unresolved".
STANDARDISE_EXACT_UNRESOLVED_TERMS = True

APPLY_SPECIES_ITALICS = False

# Y-axis is kept at 0-100 so the figure shows how much of the total profile
# is represented by the displayed top-12 feature set.
Y_AXIS_MAX = 100

# Publication sizing.
MM_TO_INCH = 1 / 25.4
FIGSIZE = (180 * MM_TO_INCH, 135 * MM_TO_INCH)
EXPORT_DPI = 600
BASE_FONT_SIZE = 8

# Colour-blind-aware palette. "Others" and "Unresolved" use fixed greys.
FEATURE_PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#D55E00",
    "#CC79A7", "#56B4E9", "#F0E442", "#332288",
    "#88CCEE", "#44AA99", "#AA4499", "#999933",
    "#117733", "#DDCC77", "#882255",
]
SPECIAL_COLOURS = {
    "Others": "#BDBDBD",
    "Unresolved": "#707070",
}

# Optional mapping of common KEGG level-3 identifiers to pathway names.
# Add or edit entries here if your table contains additional KO pathway codes.
KEGG_LEVEL3_MAP = {
    "ko00010": "Glycolysis / Gluconeogenesis",
    "ko00020": "Citrate cycle (TCA cycle)",
    "ko00190": "Oxidative phosphorylation",
    "ko00230": "Purine metabolism",
    "ko00240": "Pyrimidine metabolism",
    "ko00250": "Alanine, aspartate and glutamate metabolism",
    "ko00270": "Cysteine and methionine metabolism",
    "ko00520": "Amino sugar and nucleotide sugar metabolism",
    "ko00620": "Pyruvate metabolism",
    "ko00630": "Glyoxylate and dicarboxylate metabolism",
    "ko00640": "Propanoate metabolism",
    "ko00680": "Methane metabolism",
    "ko00720": "Carbon fixation pathways in prokaryotes",
    "ko00910": "Nitrogen metabolism",
    "ko00970": "Aminoacyl-tRNA biosynthesis",
    "ko02010": "ABC transporters",
    "ko02020": "Two-component system",
    "ko02024": "Quorum sensing",
    "ko03010": "Ribosome",
    "ko03030": "DNA replication",
    "ko03070": "Bacterial secretion system",
}


# =============================================================================
# Helper functions
# =============================================================================

def script_dir():
    """Return the directory containing this script."""
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()


def resolve_input_path():
    """Find the input CSV beside the script or in the current working directory."""
    direct = Path(INPUT_CSV)
    beside_script = script_dir() / INPUT_CSV

    if direct.exists():
        return direct
    if beside_script.exists():
        return beside_script

    raise FileNotFoundError(
        f"Could not find {INPUT_CSV}. Place the CSV in the same folder as "
        "this script or edit INPUT_CSV to the full file path."
    )


def choose_font():
    """Use Arial when available; otherwise use a safe Matplotlib fallback."""
    available = {f.name for f in font_manager.fontManager.ttflist}
    for font_name in ["Arial", "Helvetica", "DejaVu Sans"]:
        if font_name in available:
            return font_name
    return "DejaVu Sans"


def configure_matplotlib():
    """Apply journal-friendly typography and vector-font settings."""
    mpl.rcParams.update({
        "font.family": choose_font(),
        "font.size": BASE_FONT_SIZE,
        "axes.labelsize": BASE_FONT_SIZE,
        "axes.titlesize": BASE_FONT_SIZE + 1,
        "xtick.labelsize": BASE_FONT_SIZE - 1,
        "ytick.labelsize": BASE_FONT_SIZE - 1,
        "legend.fontsize": BASE_FONT_SIZE - 1,
        "legend.title_fontsize": BASE_FONT_SIZE,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.8,
    })


def remove_rank_prefixes(label):
    """
    Remove taxonomy prefixes and keep the terminal rank in lineage strings.

    Examples
    --------
    g__Acinetobacter;s__Acinetobacter wanghuae -> Acinetobacter wanghuae
    g__Unclassified;s__Actinomycetota bacterium -> Actinomycetota bacterium
    """
    if pd.isna(label):
        return ""

    text = str(label).strip()

    if text in KEGG_LEVEL3_MAP:
        return KEGG_LEVEL3_MAP[text]

    # Keep terminal feature when a full lineage is provided.
    if ";" in text or "|" in text:
        parts = re.split(r"[;|]", text)
        parts = [p.strip() for p in parts if p.strip()]
        if parts:
            text = parts[-1]

    # Remove rank prefixes at start or after delimiters/spaces.
    text = re.sub(r"(^|[\s,;|])([dkpcofgs])_{1,2}", r"\1", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def standardise_label(label):
    """Clean feature labels and standardise only true residual categories."""
    raw = "" if pd.isna(label) else str(label).strip()
    raw_low = raw.lower()

    if raw_low in {"other", "others"}:
        return "Others"

    unresolved_exact = {"unresolved", "unassigned", "unannotated", "unknown", "not annotated"}

    if STANDARDISE_EXACT_UNRESOLVED_TERMS and raw_low in unresolved_exact:
        return "Unresolved"

    if POOL_UNCLASSIFIED_AS_UNRESOLVED and "unclassified" in raw_low:
        return "Unresolved"

    clean = remove_rank_prefixes(raw)
    clean_low = clean.lower()

    if clean_low in {"other", "others"}:
        return "Others"

    if STANDARDISE_EXACT_UNRESOLVED_TERMS and clean_low in unresolved_exact:
        return "Unresolved"

    return clean


def load_table(csv_path):
    """Load a feature-by-sample abundance table."""
    df = pd.read_csv(csv_path)

    if FEATURE_COLUMN in df.columns:
        feature_col = FEATURE_COLUMN
    else:
        feature_col = df.columns[0]
        warnings.warn(
            f"Expected feature column '{FEATURE_COLUMN}' was not found. "
            f"Using first column '{feature_col}' instead."
        )

    df = df.rename(columns={feature_col: "Label"})
    sample_cols = [c for c in df.columns if c != "Label"]

    for col in sample_cols:
        df[col] = (
            df[col].astype(str)
            .str.replace(",", ".", regex=False)
            .replace({"": np.nan, "nan": np.nan, "None": np.nan})
        )
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    df["Label"] = df["Label"].apply(standardise_label)

    # Merge duplicates created by label cleaning.
    df = df.groupby("Label", as_index=False)[sample_cols].sum()

    sums = df[sample_cols].sum(axis=0)
    print("Original sample column-sum range:", sums.min(), "to", sums.max())
    if not np.allclose(sums.values, np.ones(len(sample_cols)), atol=1e-3):
        warnings.warn(
            "Columns do not sum to approximately 1. The script will not normalize "
            "them because this figure is intended to show the supplied abundance values."
        )

    return df, sample_cols


def select_top_features(df, sample_cols, top_n=12):
    """
    Select the top N features by mean abundance, without pooling omitted features.

    "Others" and "Unresolved" are included only if they are already among the
    selected top N after label cleaning. For visual consistency, they are placed
    at the top of stacked bars by plotting them last.
    """
    work = df.copy()
    work["mean_abundance"] = work[sample_cols].mean(axis=1)
    selected = work.sort_values("mean_abundance", ascending=False).head(top_n).copy()

    # Plot regular features first, then special residual categories last.
    special_order = ["Unresolved", "Others"]
    regular = selected[~selected["Label"].isin(special_order)].sort_values(
        "mean_abundance", ascending=False
    )
    special = selected[selected["Label"].isin(special_order)].copy()
    special["special_rank"] = special["Label"].map({label: i for i, label in enumerate(special_order)})
    special = special.sort_values("special_rank")

    selected = pd.concat([regular, special], ignore_index=True)
    selected = selected.drop(columns=[c for c in ["mean_abundance", "special_rank"] if c in selected.columns])

    displayed_sum = selected[sample_cols].sum(axis=0)
    coverage = pd.DataFrame({
        "Sample": sample_cols,
        "Displayed_top12_sum_fraction": displayed_sum.values,
        "Displayed_top12_sum_percent": displayed_sum.values * 100.0,
        "Omitted_fraction": 1.0 - displayed_sum.values,
        "Omitted_percent": (1.0 - displayed_sum.values) * 100.0,
    })

    return selected, coverage


def sample_compartment(sample_id):
    """Infer broad sample compartment from sample ID."""
    s = str(sample_id).upper()
    if s.startswith("GM"):
        return "GM"
    if s.startswith("RSCR"):
        return "CR"
    if s.startswith("RS"):
        return "RS"
    return "Other"


def add_group_separators(ax, sample_cols):
    """Add subtle vertical separators between GM, control rhizosphere, and RS."""
    groups = [sample_compartment(s) for s in sample_cols]
    for i in range(1, len(groups)):
        if groups[i] != groups[i - 1]:
            ax.axvline(i - 0.5, color="0.35", linewidth=0.6)


def wrap_label(label, width=32):
    """Wrap long legend labels without breaking words."""
    if label in {"Others", "Unresolved"}:
        return label
    return "\n".join(textwrap.wrap(str(label), width=width, break_long_words=False))


def is_species_like(label):
    """Heuristic for italicising species-level taxonomic names."""
    if not APPLY_SPECIES_ITALICS:
        return False

    if label in {"Others", "Unresolved"}:
        return False

    text = str(label).strip()
    text_no_cand = text.replace("Candidatus ", "")
    words = text_no_cand.split()

    if len(words) < 2:
        return False

    genus, second = words[0], words[1]

    # Avoid broad placeholder taxa.
    if second.lower() in {
        "bacterium", "archaeon", "virus", "viruses", "phage",
        "metagenome", "genome",
    }:
        return False
    if genus.lower() in {"uncultured", "unclassified", "unknown"}:
        return False
    if genus.endswith(("ota", "etes", "ia", "aceae", "ales", "ales.")):
        return False

    if not re.match(r"^[A-Z][A-Za-z-]+$", genus):
        return False

    if second == "sp.":
        return True
    if re.match(r"^[a-z][a-z-]+$", second):
        return True

    return False


def feature_colours(labels):
    """Assign colours to features."""
    regular = [x for x in labels if x not in SPECIAL_COLOURS]
    colour_map = {}
    for i, label in enumerate(regular):
        colour_map[label] = FEATURE_PALETTE[i % len(FEATURE_PALETTE)]
    for label, colour in SPECIAL_COLOURS.items():
        if label in labels:
            colour_map[label] = colour
    return [colour_map[label] for label in labels]


def save_figure(fig, output_base):
    """Save each plot as PDF, SVG, and high-resolution PNG."""
    for ext in ["pdf", "svg", "png"]:
        path = f"{output_base}.{ext}"
        if ext == "png":
            fig.savefig(path, dpi=EXPORT_DPI, bbox_inches="tight", pad_inches=0.15)
        else:
            fig.savefig(path, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def plot_bar(df_top, sample_cols, output_dir):
    """Generate the non-normalized top-12 stacked abundance plot."""
    labels = df_top["Label"].tolist()
    data = df_top.set_index("Label")[sample_cols] * 100.0
    colours = feature_colours(labels)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    fig.subplots_adjust(right=0.68)

    x = np.arange(len(sample_cols))
    bottom = np.zeros(len(sample_cols))

    for label, colour in zip(labels, colours):
        values = data.loc[label].values
        ax.bar(
            x, values, bottom=bottom, width=0.82,
            color=colour, edgecolor="white", linewidth=0.25,
            label=wrap_label(label),
        )
        bottom += values

    add_group_separators(ax, sample_cols)

    if SHOW_TITLE:
        ax.set_title(PLOT_TITLE)

    ax.set_ylabel("Relative abundance in total profile (%)")
    ax.set_xlabel("Sample")
    ax.set_xticks(x)
    ax.set_xticklabels(sample_cols, rotation=90, ha="center")
    ax.set_xlim(-0.6, len(sample_cols) - 0.4)
    ax.set_ylim(0, Y_AXIS_MAX)
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
    ax.grid(axis="y", linestyle=":", linewidth=0.45, alpha=0.6)
    ax.set_axisbelow(True)

    legend = ax.legend(
        title=LEGEND_TITLE,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        borderaxespad=0.0,
        handlelength=1.0,
        handletextpad=0.5,
    )

    for text_obj, raw_label in zip(legend.get_texts(), labels):
        if is_species_like(raw_label):
            text_obj.set_fontstyle("italic")

    save_figure(fig, output_dir / f"{FIGURE_PREFIX}_top12_non_normalized_barplot")


def main():
    configure_matplotlib()

    csv_path = resolve_input_path()
    output_dir = script_dir() / OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    df, sample_cols = load_table(csv_path)
    df_top, coverage = select_top_features(df, sample_cols, TOP_N)

    df_top.to_csv(output_dir / f"{FIGURE_PREFIX}_top12_non_normalized_plotting_table.csv", index=False)
    coverage.to_csv(output_dir / f"{FIGURE_PREFIX}_top12_displayed_coverage.csv", index=False)

    print("Input:", csv_path)
    print("Samples:", len(sample_cols))
    print("Displayed top-12 features:")
    for label in df_top["Label"].tolist():
        print(" -", label)
    print(
        "Displayed top-12 coverage range:",
        f"{coverage['Displayed_top12_sum_percent'].min():.2f}% to "
        f"{coverage['Displayed_top12_sum_percent'].max():.2f}%"
    )

    plot_bar(df_top, sample_cols, output_dir)

    print("Output folder:", output_dir)


if __name__ == "__main__":
    main()
