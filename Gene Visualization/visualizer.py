import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np
import seaborn as sns


def manhattan_plot(df, title="Manhattan Plot"):
    """
    Genom çapında ilişkilendirme çalışması (GWAS) sonuçlarını görselleştirir.

    Parametreler:
    df (pd.DataFrame): 'chromosome', 'position' ve 'p_value' sütunlarını içeren DataFrame.
    title (str): Grafiğin başlığı.
    """
    if (
        "chromosome" not in df.columns
        or "position" not in df.columns
        or "p_value" not in df.columns
    ):
        raise ValueError(
            "Manhattan Plot için DataFrame 'chromosome', 'position' ve 'p_value' sütunlarını içermelidir."
        )

    df = df.copy()
    df["neg_log_p_value"] = -1 * df["p_value"].apply(lambda x: np.log10(x))

    df["ind"] = range(len(df))
    df_grouped = df.groupby(["chromosome"])

    fig, ax = plt.subplots(figsize=(12, 6))

    colors = ["gray", "dimgray"]
    x_labels = []
    x_labels_pos = []

    for num, (name, group) in enumerate(df_grouped):
        group.plot(
            kind="scatter", x="ind", y="neg_log_p_value", color=colors[num % 2], ax=ax
        )
        x_labels.append(name)
        x_labels_pos.append((group["ind"].iloc[-1] + group["ind"].iloc[0]) / 2)

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel("Kromozom")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)
    plt.show()


def venn_diagram(sets, set_labels, title="Venn Diyagramı"):
    """
    Gen setleri arasındaki kesişimleri görselleştirir.

    Parametreler:
    sets (list of sets): Kıyaslanacak gen setlerinin bir listesi.
    set_labels (list of str): Her bir setin etiketi.
    title (str): Grafiğin başlığı.
    """
    plt.figure(figsize=(8, 8))

    if len(sets) == 2:
        venn2(subsets=sets, set_labels=set_labels)
    elif len(sets) == 3:
        venn3(subsets=sets, set_labels=set_labels)
    else:
        raise ValueError("Venn diyagramı sadece 2 veya 3 set için desteklenmektedir.")

    plt.title(title)
    plt.show()


def heatmap(df, title="Gen İfadesi Isı Haritası"):
    """
    Gen ekspresyon verilerini bir ısı haritası olarak görselleştirir.

    Parametreler:
    df (pd.DataFrame): 'Gene' (satır indeksi olarak), 'Condition' ve
                       'Expression' sütunlarını içeren DataFrame.
    title (str): Grafiğin başlığı.
    """
    if (
        "Gene" not in df.columns
        or "Condition" not in df.columns
        or "Expression" not in df.columns
    ):
        raise ValueError(
            "Heatmap için DataFrame 'Gene', 'Condition' ve 'Expression' sütunlarını içermelidir."
        )

    df_pivot = df.pivot_table(index="Gene", columns="Condition", values="Expression")

    plt.figure(figsize=(10, 8))
    sns.heatmap(df_pivot, annot=True, cmap="viridis", fmt=".2f")
    plt.title(title)
    plt.xlabel("Koşullar")
    plt.ylabel("Genler")
    plt.show()


def volcano_plot(
    df, fold_change_threshold=1.0, p_value_threshold=0.05, title="Volcano Plot"
):
    """
    Gen ekspresyon değişimlerini bir volkan grafiği olarak görselleştirir.

    Parametreler:
    df (pd.DataFrame): 'log2_fold_change' ve 'p_value' sütunlarını içeren DataFrame.
    fold_change_threshold (float): Anlamlı kabul edilecek log2 kat değişim eşik değeri.
    p_value_threshold (float): Anlamlı kabul edilecek p-değer eşik değeri.
    title (str): Grafiğin başlığı.
    """
    if "log2_fold_change" not in df.columns or "p_value" not in df.columns:
        raise ValueError(
            "Volcano Plot için DataFrame 'log2_fold_change' ve 'p_value' sütunlarını içermelidir."
        )

    df = df.copy()
    df["neg_log_p_value"] = -np.log10(df["p_value"])

    df["label"] = "Non-significant"
    df.loc[
        (df["log2_fold_change"].abs() > fold_change_threshold)
        & (df["p_value"] < p_value_threshold),
        "label",
    ] = "Significant"

    colors = {"Non-significant": "grey", "Significant": "red"}

    plt.figure(figsize=(10, 8))
    plt.scatter(
        x=df["log2_fold_change"],
        y=df["neg_log_p_value"],
        c=df["label"].map(colors),
        alpha=0.6,
        s=10,
    )

    plt.axvline(x=fold_change_threshold, color="black", linestyle="--", linewidth=1)
    plt.axvline(x=-fold_change_threshold, color="black", linestyle="--", linewidth=1)
    plt.axhline(
        y=-np.log10(p_value_threshold), color="black", linestyle="--", linewidth=1
    )

    plt.title(title)
    plt.xlabel("log2 Kat Değişimi")
    plt.ylabel("-log10(p-value)")
    plt.show()
