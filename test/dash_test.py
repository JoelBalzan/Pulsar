from dash import Dash, dcc, callback, Input, Output
import plotly.express as px
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import numpy as np

df = np.load('/media/joel/JB_EXHDD/J1809-1943/P970/P970_Spectra_v2_704_4032_on_2.npy')

app = Dash(__name__)

app.layout = dbc.App([
    dbc.Header("Dynamic Spectrum"),

    dbc.Card(
        dbc.Graph(figure=px.imshow(a, aspect='auto', origin='lower'),width=50)
    )
])