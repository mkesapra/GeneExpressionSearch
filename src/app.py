#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 10:52:38 2024

@author: kesaprm
"""
import pathlib
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import pandas as pd
import numpy as np


def load_data(data_file: str) -> pd.DataFrame:
    '''
    Load data from /data directory
    '''
    PATH = pathlib.Path(__file__).parent
    DATA_PATH = PATH.joinpath("data").resolve()
    return pd.read_excel(DATA_PATH.joinpath(data_file))

#df = pd.read_excel("SomeDaVinciGenes.xlsx")
df = load_data("SomeDaVinciGenes.xlsx")

app = dash.Dash(__name__)
# Declare server for Heroku deployment. Needed for Procfile.
server = app.server


# Dropdown options for gene names
gene_options = [{'label': gene, 'value': gene} for gene in df.columns[1:-4]]

app.layout = html.Div([
    html.H1("Gene Expression Lookup", style={'text-align': 'center', 'margin-bottom': '20px'}),  # Title at the top
    html.Div([
        dcc.Dropdown(
            id='gene-dropdown',
            options=gene_options,
            placeholder="Select a gene",
            clearable=True,
            style={'width': '300px'}  # Set a fixed width for the dropdown
        ),
        html.Button('Look up', id='lookup-button', style={'margin-left': '5px', 'height': '35px'})  # Add margin to button
    ], style={'display': 'flex', 'align-items': 'center'}),  # Flexbox for alignment
    dcc.RadioItems(
        id='data-type',
        options=[
            {'label': 'Raw Data', 'value': 'raw'},
            {'label': 'Log2(1 + expression)', 'value': 'log'}
        ],
        value='raw',
        labelStyle={'display': 'inline-block'},
        style={'margin-top': '15px'},
    ),
    dcc.Graph(id='gene-plot')
])

@app.callback(
    Output('gene-plot', 'figure'),
    [Input('lookup-button', 'n_clicks'),
     Input('data-type', 'value')],
    [State('gene-dropdown', 'value')]
)
def update_plot(n_clicks, data_type, gene_name):
    if n_clicks is None or not gene_name:
        raise dash.exceptions.PreventUpdate  # Prevent any plot from loading initially

    gene_columns = {col.lower(): col for col in df.columns[1:-4]}  # Map lowercase to original columns

    gene_name_lower = gene_name.lower()
    if gene_name_lower in gene_columns:
        gene_col = gene_columns[gene_name_lower]
        
        # Select relevant data based on 'Edge/Center' and 'gene_col'
        df_single_gene = df[['Day', 'Edge/Center', gene_col]]
        df_single_gene = df_single_gene[df_single_gene['Edge/Center'].isin(['c', 'e'])]

        # Conditionally apply log transformation and calculate max y-axis limit for log data
        if data_type == 'log':
            df_single_gene[gene_col] = np.log2(1 + df_single_gene[gene_col])
            max_expression_value = np.log2(1 + df.iloc[:, 1:-4].max().max())
        else:
            max_expression_value = None  # No y-axis limit for raw data

        # Reshape data and calculate mean lines
        df_single_gene = df_single_gene.melt(id_vars=['Day', 'Edge/Center'],
                                             value_vars=gene_col,
                                             var_name='Gene',
                                             value_name='Expression')

        mean_edge = df_single_gene[df_single_gene['Edge/Center'] == 'e'].groupby('Day')['Expression'].mean().reset_index()
        mean_center = df_single_gene[df_single_gene['Edge/Center'] == 'c'].groupby('Day')['Expression'].mean().reset_index()

        # Create scatter plot with mean lines
        fig = go.Figure()

        for label, color, label_name in [('e', 'red', 'Edge'), ('c', 'blue', 'Center')]:
            filtered_data = df_single_gene[df_single_gene['Edge/Center'] == label]
            fig.add_trace(go.Scatter(
                x=filtered_data['Day'], 
                y=filtered_data['Expression'],
                mode='markers',
                marker=dict(color=color),
                name=f'{label_name}'
            ))

        fig.add_trace(go.Scatter(x=mean_edge['Day'], y=mean_edge['Expression'],
                                 mode='lines+markers',
                                 line=dict(color='red', dash='dash'),
                                 name='Mean Edge Expression'))

        fig.add_trace(go.Scatter(x=mean_center['Day'], y=mean_center['Expression'],
                                 mode='lines+markers',
                                 line=dict(color='blue', dash='dash'),
                                 name='Mean Center Expression'))

        # Apply y-axis limit only for log-transformed data
        yaxis_range = [0, max_expression_value] if data_type == 'log' else None

        fig.update_layout(title=f'Expression of {gene_col} Over Days ({data_type.title()})',
                          xaxis_title='Day',
                          yaxis_title='Expression' if data_type == 'raw' else 'Log2(1 + Expression)',
                          xaxis=dict(tickvals=[0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 16, 19, 21]),
                          yaxis=dict(range=yaxis_range),  # Set y-axis range conditionally
                          height=600)

        return fig
    else:
        return go.Figure()  # Return an empty figure if gene name is invalid


if __name__ == '__main__':
    app.run_server(debug=False)
