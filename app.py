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

st.title("Annexin Receptor Dashboard")

st.sidebar.header("Gene Search")

gene = st.sidebar.text_input(
    "Enter Gene Symbol",
    value="ANXA1"
)

st.sidebar.markdown("---")

st.sidebar.write(
    f"Cells: {adata.n_obs:,}"
)

st.sidebar.write(
    f"Genes: {adata.n_vars:,}"
)

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
                Mean_Expression=("expression","mean"),
                Percent_Positive=("expression",
                                  lambda x:(x>0).mean()*100)
            )
            .sort_values(
                "Mean_Expression",
                ascending=False
            )
        )

        st.subheader(
            "Top 10 Expressing Cell Types"
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

        st.subheader(
            "Expression Across Cell Types"
        )

        fig, ax = plt.subplots(
            figsize=(10,5)
        )

        summary["Mean_Expression"].plot.bar(
            ax=ax
        )

        plt.tight_layout()

        st.pyplot(fig)

        st.subheader(
            "Condition Comparison"
        )

        cond_df = pd.DataFrame({
            "Condition": adata.obs["Condition"],
            "expression": expr
        })

        cond_summary = (
            cond_df.groupby("Condition")
            ["expression"]
            .mean()
        )

        fig2, ax2 = plt.subplots()

        cond_summary.plot.bar(
            ax=ax2
        )

        plt.tight_layout()

        st.pyplot(fig2)
