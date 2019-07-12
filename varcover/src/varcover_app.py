import os
import sys
import random
import gzip
import base64
import datetime
import urllib
import subprocess
import io
from io import StringIO
from pathlib import Path

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_table as dt

import pandas as pd
import plotly.graph_objs as go


home = str(Path.home())
sys.path.append(home + '/docker/varcover/varcover/src/')
from varcover_ensembl import *
from varcover_preprocess import *
from varcover import *
from varcover_app_utils import parse_contents
from varcover_app_utils import parse_vcf

sys.path.append(home + '/git/pandasVCF')
from pandasvcf import *

from flask import Flask

#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
#external_stylesheets = ['./data/bWLwgP.css']
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css', 'https://codepen.io/chriddyp/pen/brPBPO.css']


server = Flask(__name__)

app = dash.Dash(__name__,server=server, external_stylesheets=external_stylesheets)
#app.config['suppress_callback_exceptions'] = True

colors = {
          'background': '#111111',
          'text': '#7FDBFF'
          }


app.layout = html.Div(children=[

  ### HEADER
    html.Div(children = [
            html.H1(children='Icahn School of Medicine at Mount Sinai',
            style={'background-image': 'linear-gradient(to right,#00aeef,#d80b8c)',
                  'color':'white',
                  'font-size':22,
                  'marginBottom': '0.0em',
                  'marginTop': '1.0em',
                  'padding-left':'5px'}),

            html.H1(children='Department of Genetics & Genomic Sciences',
            style={'background-image': 'linear-gradient(to right,#00aeef,#d80b8c)',
                  'color':'white',
                  'font-size':22,
                  'marginBottom': '4.0em',
                  'marginTop': '0.0em',
                  'textAlign':'right',
                  'padding-right':'5px'})
                  ]),



  ### WELCOME & CAVEATS
    html.Div(children = [
            html.H1(children='Welcome to the VarCover Web App!',
            style={'textAlign': 'center',
                   'color': '#d80b8c',
                  'marginBottom': '1.5em'}),

    html.H5("Important Note: This is a public demonstration site and provides "
            "No Privacy Or Security Features.",
            style={'textAlign':'center', 'color':'#00aeef',
                  'marginBottom': '0.0em'}),

    html.H6("Consider using the VarCover package for tasks that require greater "
            "data privacy, security, or compute resources: https://github.com/erscott/varcover",
            style={'textAlign':'center', 'color':'#00aeef',
                  'marginTop': '0.0em',
                  'marginBottom': '5.0em'})
                  ]),



  ### INTERACTIVE COMPONENTS

    html.H2('Please select:',
            style={'textAlign':'center', 'color':'black',
                  'padding-left':'0px'}),

  ### DROP DOWN MENU COST METRIC
    html.Div(children = [
        html.H3(' Weighting Metric', style={'textAlign': 'center',
                'color':'#d80b8c', 'marginBottom': '0.0em'}),

        html.Abbr("hover for help", title="Allele frequency logit "
                  "weighting increases the likelihood of additional alleles "
                  "in the solution set with ~ similar sample sizes.")
                  ],
                  style={'textAlign':'center', 'color':'blue', 'fontsize':40,
                       'marginTop': '1.0em'}),

        dcc.RadioItems(
            id='cost-metric',
            options=[
                {'label': 'Standard', 'value': 'standard'},
                {'label': 'Allele Frequency Logit', 'value': 'logit'},
                ],
             value='standard',
             style={'textAlign':'center',  'font-size':20,
                    'marginBottom': '2.0em'})
           ,


  ### DROP DOWN MENU REDUCE SINGLETONS

    html.Div(children = [
        html.H3('Reduce Singletons Speedup', style={'textAlign': 'center',
                    'color':'#d80b8c', 'marginBottom': '0.0em'}),

        html.Abbr("hover for help", title="Reduces computational "
                    "complexity by selecting all samples with singleton alleles "
                    "prior to solving the min-set cover problem. Most helpful "
                    "for target sets with >200 variants.")
                    ],
                    style={'textAlign':'center', 'color':'blue', 'fontsize':40,
                           'marginTop': '1.0em'}),

        dcc.RadioItems(
            id='reduce-singletons',
            options=[
                {'label': 'True', 'value': 'True'},
                {'label': 'False', 'value': 'False'}],
            value='False',
            style={'textAlign':'center',  'font-size':20,
                'marginBottom':'2.5em'}),



  ###  SUBMIT RSID OR VCF
    html.H2('Finally, submit an RSID or VCF file for analysis:',
            style={'textAlign':'center', 'color':'black',
                  'padding-left':'0px'}),


    html.Div(children = [

            html.H3('RSID file for 1KG Phase 3 Analysis',
            style={'textAlign': 'center','marginBottom': '0.0em',
                   'color': '#d80b8c'}),
            html.Abbr("file format ex.", title="#1 rsid per line: .tsv, .csv, "
                                           ".txt file\nrs6500\nrs2032582\nrs56116432")
            ],

            style={'textAlign':'center', 'color':'blue', 'fontsize':40,
                               'marginTop': '0.0em'}),

  ### UPLOAD RSID FILE
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'RSID File: ',
                'Drag and Drop or ',
                html.A('Select Files')
                ]),
             style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '0px'},
         # Allow multiple files to be uploaded
         multiple=False
         ),


    ### Displaying the uploaded results
    html.Div(id='output-data-upload'),


  ###  RSID OR VCF

    html.H3('---',
            style={'textAlign':'center', 'color':'black'}),

    html.H3('OR',
            style={'textAlign': 'center',
                   'color': '#d80b8c'}),
    html.H3('---',
            style={'textAlign':'center', 'color':'black'}),

    html.H3('Submit a VCF file',
            style={'textAlign': 'center',
                   'color': '#d80b8c'}),

    html.H6("VCF files: Please do not submit files with private health information "
            "or other privacy restrictions",
            style={'textAlign': 'center',
                   'color': 'black'}),


  ### UPLOAD VCF FILE
    dcc.Upload(
            id='upload-vcf',
            children=html.Div([
                'VCF(.gz) file: ',
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '0px'
            },
            # Allow multiple files to be uploaded
            multiple=False,
            max_size=1000000
        ),

    ### Displaying the uploaded results

    html.Div(id='output-vcf-upload'),



  ### ATTRIBUTIONS AND FOOTER

    html.Div(children = [
        html.Div(html.H3(
                dcc.Link('This website hosts 1000 Genomes Project Phase 3 data',
                         href=f'http://www.internationalgenome.org/'),
            style={'textAlign':'center', 'color':'black',
                  'marginBottom':'4.0em',
                  'marginTop':'3.0em'})),

        html.H1(children='Icahn School of Medicine at Mount Sinai',
            style={'background-image': 'linear-gradient(to right,#00aeef,#d80b8c)',
                  'color':'white',
                  'font-size':22,
                  'marginBottom': '0.0em',
                  'marginTop': '2.0em',
                  'padding-left':'5px'}),

        html.H1(children='Department of Genetics & Genomic Sciences',
            style={'background-image': 'linear-gradient(to right,#00aeef,#d80b8c)',
                  'color':'white',
                  'font-size':22,
                  'marginBottom': '0.0em',
                  'marginTop': '0.0em',
                  'textAlign':'right',
                  'padding-right':'5px'}),
        ])

])



### APP CALLBACKS

# handle uploaded rsid file and pass to table parser
@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified'),
               State('cost-metric','value'),
               State('reduce-singletons','value')
              ])
def update_output(list_of_contents, list_of_names, list_of_dates,
                 cost_metric, reduce_singletons):
    if list_of_contents is not None:
        rsid_html = parse_contents(list_of_contents, list_of_names, list_of_dates,
                                   cost_metric, reduce_singletons)
        print('Returning children, type={}'.format(type(rsid_html)))
        return rsid_html

    else:
        rsid_html = html.Div([
            html.H5('Awaiting RSID Upload',
                    style={'textAlign': 'center',
                           'color': 'blue'}),
            html.H6('(Max Number of RSID is 500)',
                    style={'textAlign': 'center',
                           'color': 'blue'})])

        return rsid_html


# handle uploaded vcf file and pass to table parser
@app.callback(Output('output-vcf-upload', 'children'),
              [Input('upload-vcf', 'contents')],
              [State('upload-vcf', 'filename'),
               State('upload-vcf', 'last_modified'),
               State('cost-metric','value'),
               State('reduce-singletons','value')
              ])
def update_output(list_of_contents, list_of_names, list_of_dates,
                 cost_metric, reduce_singletons):
    if list_of_contents is not None:
        vcf_html = parse_vcf(list_of_contents, list_of_names, list_of_dates,
                            cost_metric, reduce_singletons)

        return vcf_html

    else:
        vcf_html = html.Div([
            html.H5('Awaiting VCF Upload',
                    style={'textAlign': 'center',
                           'color': 'blue'}),
            html.H6('(Max File Size is 1MB)',
                    style={'textAlign': 'center',
                           'color': 'blue'})])

        return vcf_html


if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=5000)
