# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.


import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output

app = Dash(__name__)

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options



app = Dash(__name__)


app.layout = html.Div([
    html.H4("Analysis of the restaurant's revenue"),
    html.P("x-axis:"),
    dcc.Checklist(
        id='x-axis', 
        options=['status'],
        value=['status'], 
        inline=True
    ),
    html.P("y-axis:"),
    dcc.RadioItems(
        id='y-axis', 
        options=['ENST00000461347.5','ENST00000322723.8','ENST00000484328.1','ENST00000454824.5','ENST00000417652.5','ENST00000453992.5','ENST00000436894.1','ENST00000494618.1','ENST00000356936.5','ENST00000466274.1'],
        value='ENST00000461347.5', 
        inline=True
    ),
    dcc.Graph(id="graph"),
])


@app.callback(
    Output("graph", "figure"), 
    Input("x-axis", "value"), 
    Input("y-axis", "value"))
def generate_chart(x, y):
    INPUT_PATH= "/data/projects/shared_data/collab_data/subsample/boxplot_data"
    INPUT_FILENAME = "lgg_gbm_gtex.csv"
    df = pd.read_csv(INPUT_PATH+"/"+INPUT_FILENAME, sep= ",", decimal=',', index_col=0, header =0)
    df = df.T
    df = df.astype(float)
    df_pheno = pd.read_csv("/data/projects/shared_data/collab_data/subsample/boxplot_data/lgg_gbm_gtex_pheno.csv", sep= ",")
    df['status'] = df_pheno['status'].to_list()
    fig = px.box(df, x=x, y=y)
    fig.update_yaxes(nticks=10, tickformat='.%2f')
    return fig


app.run_server(debug=True, port=1269)


# fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")

# app.layout = html.Div(children=[
#     html.H1(children='Hello Dash'),

#     html.Div(children='''
#         Dash: A web application framework for your data.
#     '''),

#     dcc.Graph(
#         id='example-graph',
#         figure=fig
#     )
# ])

# if __name__ == '__main__':
#     app.run_server(debug=True, port=8052)