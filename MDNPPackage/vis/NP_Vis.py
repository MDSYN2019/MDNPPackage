import pandas as pd
import numpy as np
import sys

sys.path.append("/home/sang/Desktop/git/MDNPPackage/MDNPPackage")

# plotly functionalities 
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots
import core.NP_Core as core 
import streamlit as st

def dashboard_NP(dataframe: pd.DataFrame):
    """
    first function for dashboard
    """
    fig = px.scatter_3d(dataframe, x='X', y='Y', z='Z', color='NAME')
    return fig 

def return_streamlit_NP(dataframe: pd.DataFrame):
    plotly_NP = dashboard_NP(dataframe)
    st.plotly_chart(plotly_NP, use_container_width=True)

