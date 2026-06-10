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

# ==================================================
# SIDEBAR
# ==================================================

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

st.sidebar.write(
    f"Cells: {adata.n_obs:,}"
)

st.sidebar.write(
    f"Genes: {adata.n_vars:,}"
)

st.sidebar.markdown("---")

cell_types = sorted(
    adata.obs["cell_type"].unique()
)

selected_cell_types = st.sidebar.multiselect(
    "Cell Type Filter",
    cell_types,
    default=cell_types
)

conditions = sorted(
    adata.obs["Condition"].unique()
)

selected_conditions = st.sidebar.multiselect(
    "Condition Filter",
    conditions,
    default=conditions
)

# ==================================================
# GENE ANALYSIS
# ==================================================

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
            "Condition": adata.obs["Condition"],
            "expression": expr
        })

        df = df[
            df["cell_type"].isin(selected_cell_types)
        ]

        df = df[
            df["Condition"].isin(selected_conditions)
        ]

        if len(df) == 0:

            st.warning(
                "No cells remain after filtering."
            )

            st.stop()

        cond_df = df.copy()

        summary = (
            df.groupby("cell_type")
             .agg(
                Mean_Expression=("expression", "mean"),

                Percent_Positive=(
                    "expression",
                     lambda x: (x > 0).mean() * 100
                ),

                Cell_Count=(
                    "expression",
                    "count"
                )
            )
            .sort_values(
                "Mean_Expression",
                ascending=False
            )
        )

        if summary.empty:

            st.warning(
                "No cells remain after the selected filters."
            )

            st.stop()

        highest_cell_type = summary.index[0]

        highest_condition = (
            cond_df.groupby("Condition")
            ["expression"]
            .mean()
            .idxmax()
        )

        mean_expression = df["expression"].mean()

        positive_pct = (
            (df["expression"] > 0).mean()
            * 100
        )
        # ==========================================
        # GENE SUMMARY
        # ==========================================

        st.subheader("Gene Summary")

        col1, col2, col3, col4 = st.columns(4)

        col1.metric(
            "Mean Expression",
            f"{mean_expression:.3f}"
        )

        col2.metric(
            "% Positive Cells",
            f"{positive_pct:.1f}%"
        )

        col3.metric(
            "Highest Cell Type",
            highest_cell_type
        )

        col4.metric(
            "Highest Condition",
            highest_condition
        )

        # ==========================================
        # LAB SUTTI INSIGHT
        # ==========================================

        st.success(
            f"""
        🧠 Lab Sutti Insight

        {gene} shows strongest enrichment in {highest_cell_type}.

        Approximately {positive_pct:.1f}% of analysed cells express this gene.

        The highest condition-level expression is observed in {highest_condition}.

        This suggests that {highest_cell_type} may be a major contributor to the biological signal associated with {gene}.
        """
        )
        
        # ==========================================
        # TOP CELL TYPES TABLE
        # ==========================================

        st.subheader("Top 10 Expressing Cell Types")

        st.info(
            """
Mean Expression = average expression across all cells within a cell type.

Percent Positive = percentage of cells expressing the gene (>0 expression).

Cell types are ranked by mean expression.
"""
        )

        st.dataframe(
            summary.head(10)
        )

        csv = summary.to_csv()

        st.download_button(
            "Download CSV",
            csv,
            file_name=f"{gene}_summary.csv"
        )

        # ==========================================
        # CELL TYPE BARPLOT
        # ==========================================

        st.subheader(
            "Expression Across Cell Types"
        )

        fig, ax = plt.subplots(
            figsize=(10, 5)
        )

        summary["Mean_Expression"].plot.bar(
            ax=ax
        )

        ax.set_ylabel(
            "Mean Expression"
        )

        ax.set_xlabel(
            "Cell Type"
        )

        plt.tight_layout()

        st.pyplot(fig)

        # ==========================================
        # CONDITION COMPARISON
        # ==========================================

        st.subheader(
            "Condition Comparison"
        )

        cond_summary = (
            cond_df.groupby("Condition")
            ["expression"]
            .mean()
        )

        fig2, ax2 = plt.subplots(
            figsize=(6, 5)
        )

        cond_summary.plot.bar(
            ax=ax2
        )

        ax2.set_ylabel(
            "Mean Expression"
        )

        ax2.set_xlabel(
            "Condition"
        )

        plt.tight_layout()

        st.pyplot(fig2)

        # ==========================================
        # DOWNLOAD FILTERED DATA
        # ==========================================

        st.subheader(
            "Download Filtered Dataset"
        )

        filtered_csv = df.to_csv(
            index=False
        )

        st.download_button(
            "Download Filtered Data",
            filtered_csv,
            file_name=f"{gene}_filtered_data.csv"
        )

        # ==========================================
        # UMAP SECTION
        # ==========================================

        st.header(
            "UMAP Visualization"
        )

        if umap_mode == "Annotated Cell Types":

            fig3, ax3 = plt.subplots(
                figsize=(10, 8)
            )

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

            displayed_cells = adata.n_obs

        else:

            temp = adata.copy()

            expr_temp = temp[:, gene].X

            if not isinstance(
                expr_temp,
                np.ndarray
            ):
                expr_temp = expr_temp.toarray()

            expr_temp = expr_temp.flatten()

            if umap_mode == "All Cells":

                st.write(
                    """
Showing all cells coloured according to gene expression intensity.
"""
                )

            elif umap_mode == "Positive Cells Only":

                temp = temp[
                    expr_temp > 0
                ]

                st.write(
                    """
Showing only cells with detectable expression (>0).
"""
                )

            elif umap_mode == "Top 10% Expressing Cells":

                cutoff = np.percentile(
                    expr_temp,
                    90
                )

                temp = temp[
                    expr_temp >= cutoff
                ]

                st.write(
                    """
Showing only the highest 10% expressing cells.
"""
                )

            elif umap_mode == "Top 50% Expressing Cells":

                cutoff = np.percentile(
                    expr_temp,
                    50
                )

                temp = temp[
                    expr_temp >= cutoff
                ]

                st.write(
                    """
Showing only the highest 50% expressing cells.
"""
                )

            fig3, ax3 = plt.subplots(
                figsize=(10, 8)
            )

            sc.pl.umap(
                temp,
                color=gene,
                frameon=False,
                show=False,
                ax=ax3
            )

            displayed_cells = temp.n_obs

        colA, colB = st.columns(2)

        colA.metric(
            "Cells Displayed",
            f"{displayed_cells:,}"
        )

        colB.metric(
            "% Dataset Displayed",
            f"{100 * displayed_cells / adata.n_obs:.1f}%"
        )

        st.pyplot(fig3)

        fig3.savefig(
            "UMAP.png",
            dpi=300,
            bbox_inches="tight"
        )

        with open(
            "UMAP.png",
            "rb"
        ) as f:

            st.download_button(
                "Download UMAP",
                f,
                file_name=f"{gene}_UMAP.png"
            )
