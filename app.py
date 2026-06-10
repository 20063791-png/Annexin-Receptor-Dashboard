
import streamlit as st

st.title("Annexin Receptor Dashboard")

st.write("Dashboard is working!")

gene = st.text_input("Enter Gene Symbol")

if gene:
    st.write("You searched:", gene)
