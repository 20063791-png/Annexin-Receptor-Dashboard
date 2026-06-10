import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gdown

st.set_page_config(
    page_title="Annexin Receptor Dashboard",
    layout="wide"
)

@st.cache_resource
def load_data():

    file_id = "1uOV7cSQr3BUHILVbC3_x4Y1NvUGFs04a"

    gdown.download(
        f"https://drive.google.com/uc?id={file_id}",
        "dataset.h5ad",
        quiet=True
    )

    adata = sc.read_h5ad("dataset.h5ad")

    return adata

adata = load_data()

# --------------------------------------------------
# SIDEBAR
# --------------------------------------------------

st.title("Annexin Receptor Dashboard")

st.sidebar.header("Gene Search")

gene = st.sidebar.text_input(
    "Enter Gene Symbol",
    value="ANXA1"
)

st.sidebar.markdown("---")

umap_mode = st.sidebar.radio(
    "UMAP Display Mode",
    [
        "Annotated Cell Types",
        "All Cells",
        "Positive Cells Only",
        "Top 10% Expressing Cells",
        "Top 50% Expressing Cells"
    ]
)

st.sidebar.markdown("---")

st.sidebar.write(f"Cells: {adata.n_obs:,}")
st.sidebar.write(f"Genes: {adata.n_vars:,}")

# --------------------------------------------------
# GENE ANALYSIS
# --------------------------------------------------

if gene:

    if gene not in adata.var_names:

        st.error(f"{gene} not found")

    else:

        st.header(f"Gene: {gene}")

        expr = adata[:, gene].X

        if not isinstance(expr, np.ndarray):
            expr = expr.toarray()

        expr = expr.flatten()

        df = pd.DataFrame({
            "cell_type": adata.obs["cell_type"],
            "expression": expr
        })

        summary = (
            df.groupby("cell_type")
            .agg(
                Mean_Expression=("expression", "mean"),
                Percent_Positive=(
                    "expression",
                    lambda x: (x > 0).mean() * 100
                )
            )
            .sort_values(
                "Mean_Expression",
                ascending=False
            )
        )

        st.subheader("Top 10 Expressing Cell Types")

        st.info(
            """
Mean Expression = average expression across all cells within a cell type.

Percent Positive = percentage of cells expressing the gene (>0 expression).

Cell types are ranked by mean expression.
"""
        )

        st.dataframe(summary.head(10))

        csv = summary.to_csv()

        st.download_button(
            "Download CSV",
            csv,
            file_name=f"{gene}_summary.csv"
        )

        # ------------------------------------------
        # CELL TYPE BARPLOT
        # ------------------------------------------

        st.subheader("Expression Across Cell Types")

        fig, ax = plt.subplots(figsize=(10, 5))

        summary["Mean_Expression"].plot.bar(ax=ax)

        ax.set_ylabel("Mean Expression")
        ax.set_xlabel("Cell Type")

        plt.tight_layout()

        st.pyplot(fig)

        # ------------------------------------------
        # CONDITION COMPARISON
        # ------------------------------------------

        st.subheader("Condition Comparison")

        cond_df = pd.DataFrame({
            "Condition": adata.obs["Condition"],
            "expression": expr
        })

        cond_summary = (
            cond_df.groupby("Condition")
            ["expression"]
            .mean()
        )

        fig2, ax2 = plt.subplots(figsize=(6, 5))

        cond_summary.plot.bar(ax=ax2)

        ax2.set_ylabel("Mean Expression")
        ax2.set_xlabel("Condition")

        plt.tight_layout()

        st.pyplot(fig2)

        # ------------------------------------------
        # UMAP
        # ------------------------------------------

        st.header("UMAP Visualization")

        if umap_mode == "Annotated Cell Types":

            fig3, ax3 = plt.subplots(figsize=(10, 8))

            sc.pl.umap(
                adata,
                color="cell_type",
                legend_loc="on data",
                frameon=False,
                show=False,
                ax=ax3
            )

            st.write(
                """
Annotated UMAP showing all cell populations.

Labels correspond to cell-type annotations.
"""
            )

        else:

            temp = adata.copy()

            expr_temp = temp[:, gene].X

            if not isinstance(expr_temp, np.ndarray):
                expr_temp = expr_temp.toarray()

            expr_temp = expr_temp.flatten()

            if umap_mode == "All Cells":

                st.write(
                    """
Showing all cells coloured according to gene expression intensity.
"""
                )

            elif umap_mode == "Positive Cells Only":

                temp = temp[expr_temp > 0]

                st.write(
                    """
Showing only cells with detectable expression (>0).
"""
                )

            elif umap_mode == "Top 10% Expressing Cells":

                cutoff = np.percentile(expr_temp, 90)

                temp = temp[expr_temp >= cutoff]

                st.write(
                    """
Showing only the highest 10% expressing cells.

This threshold is calculated across all cells in the dataset.
"""
                )

            elif umap_mode == "Top 50% Expressing Cells":

                cutoff = np.percentile(expr_temp, 50)

                temp = temp[expr_temp >= cutoff]

                st.write(
                    """
Showing only the highest 50% expressing cells.

This threshold is calculated across all cells in the dataset.
"""
                )

            fig3, ax3 = plt.subplots(figsize=(10, 8))

            sc.pl.umap(
                temp,
                color=gene,
                frameon=False,
                show=False,
                ax=ax3
            )

        st.pyplot(fig3)

        fig3.savefig(
            "UMAP.png",
            dpi=300,
            bbox_inches="tight"
        )

        with open("UMAP.png", "rb") as f:

            st.download_button(
                "Download UMAP",
                f,
                file_name=f"{gene}_UMAP.png"
            )
