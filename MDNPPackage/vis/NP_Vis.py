import pandas as pd
import numpy as np
import sys

sys.path.append("/home/sang/Desktop/git/MDNPPackage/MDNPPackage")

# plotly functionalities 
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots
import core.NP_Core as core 
#from core import CentralCoreGenerator 
# streamlit dashboard
import streamlit as st

gro_path = "/home/sang/Desktop/git/Martini3-small-molecules/models/gros" 
ligand_1_path = gro_path + "/" + "1MIMI.gro" # First ligand 
ligand_2_path = gro_path + "/" + "2NITL.gro" # Second ligand 
itp_1_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/1MIMI_cog.itp'
itp_2_path = '/home/sang/Desktop/git/Martini3-small-molecules/models/itps/cog-mono/2NITL_cog.itp'

CoreGen = core.CentralCoreGenerator(10, 50, [ligand_1_path, ligand_2_path], ['N1', 'N1'] , ['R3', 'R4'],  itp_1_path, itp_2_path, option = 'Janus')
DF = CoreGen.return_ordered_coordinates()

def dashboard_NP(dataframe: pd.DataFrame):
    """
    first function for dashboard
    """
    fig = px.scatter_3d(dataframe, x='X', y='Y', z='Z', color='NAME')
    return fig 

def return_streamlit_NP(dataframe: pd.DataFrame):
    plotly_NP = dashboard_NP(dataframe)
    st.plotly_chart(plotly_NP, use_container_width=True)

