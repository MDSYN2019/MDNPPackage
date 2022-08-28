import pandas as pd
import numpy as np
import sys

# plotly functionalities 
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

# streamlit dashboard
import streamlit as st


def dashboard_NP(dataframe: pd.DataFrame):
    """
    first function for dashboard
    """
    fig = px.scatter_3d(dataframe, x='X', y='Y', z='Z', color='name')
    return fig 
    


