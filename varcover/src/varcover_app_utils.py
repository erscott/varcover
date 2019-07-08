
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

import pandas as pd

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_table_experiments as dt

import plotly.graph_objs as go


home = str(Path.home())
sys.path.append(home + '/docker/varcover/varcover/src/')
from varcover_ensembl import *
from varcover_preprocess import *
from varcover import *

sys.path.append(home + '/git/pandasVCF')
from pandasvcf import *

from flask import Flask


### FUNCTIONS FOR UPLOADING AND TABLE RENDERING

# parse uploaded file into table
def parse_contents(contents, filename, date,
                  cost_metric, reduce_singletons):

    if reduce_singletons == 'True':
        reduce_singletons = True
    else:
        reduce_singletons = False

    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if '.csv' in filename:
            # Assume that the user uploaded a CSV file
            rsids = pd.read_csv(
                             io.StringIO(decoded.decode('utf-8')),
                             header=None,
                             comment='#')
            rsids.columns = ['rsid']
        elif '.tsv' in filename:
            # Assume that the user uploaded a CSV file
            rsids = pd.read_csv(
                               io.StringIO(decoded.decode('utf-8')),
                               sep='\t',
                               header=None,
                               comment='#')
            rsids.columns = ['rsid']
        elif '.txt' in filename:
            # Assume that the user uploaded a CSV file
            rsids = pd.read_csv(
                               io.StringIO(decoded.decode('utf-8')),
                               sep='\t',
                               header=None,
                               comment='#')
            rsids.columns = ['rsid']
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            rsids = pd.read_excel(io.BytesIO(decoded))
            rsids.columns = ['rsid']

    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    ensembl = EnsemblVar(rsids['rsid'].tolist(), genotypes=True)
    gts = ensembl.get_set_cover_df()

    vc = varcover(gts)
    vc.getCoverSet(cost=cost_metric,
                   reduceSingletons=reduce_singletons)

    vc_soln = vc.solution.reset_index()

    df = vc_soln.sort_values(['CHROM', 'POS'])

    missing_rsids = list(set(rsids['rsid'].tolist()) - set(vc_soln.query_rsid.values))

    solution_samples = vc.sample_target_allele_cnt.reset_index()


    div = [
            html.H3('VarCover Results for: {}'.format(filename),
                    style={'textAlign': 'left',
                           'color': '#00aeef'}),

            html.H6('Weighting Metric: {}'.format(cost_metric),
                    style={'textAlign': 'left',
                           'color': 'black'}),

            html.H6('Reduce Singletons: {}'.format(reduce_singletons),
                    style={'textAlign': 'left',
                           'color': 'black'}),

            html.H6(datetime.datetime.fromtimestamp(date)),
            html.H4('VarCover Solution Samples',
                    style={'textAlign': 'left',
                           'color': '#d80b8c',
                           'marginBottom':'0.0em'}),
            html.H4('n={}'.format(solution_samples.shape[0]),
                    style={'textAlign': 'left',
                           'color': '#d80b8c',
                           'marginTop':'0.0em'}),
            dt.DataTable(rows=solution_samples.to_dict('records'),
                         filterable=True),

            html.H4('VarCover Solution Matrix',
                    style={'textAlign': 'left',
                       'color': '#d80b8c'}),
            dt.DataTable(rows=df.to_dict('records'), filterable=True)
           ]

    if len(missing_rsids) == 0:

        div.extend([html.H3('No Uncovered Variants'),
                    html.Br()
                   ])

    else:

        missing_rsids = pd.DataFrame(missing_rsids, columns=['rsid'])


        div.extend([html.H4('Uncovered Variants',
                            style={'textAlign': 'left',
                                        'color': '#d80b8c'}),
                    dt.DataTable(rows=missing_rsids.to_dict('records'),
                                       filterable=True),
                    html.Br()
                   ])


    div.extend([html.A(
                    html.Button('Download RSID Solution Samples',
                                className='container',
                                style={'color':'blue',
                                       'textAlign':'center',
                                       'width':'25%'}),
                        id='download-rsid-sample-link',
                        download="VarCover_RSID_sample_set.tsv",
                        href="",
                        target="_blank"),

                html.Br(),

                html.A(
                     html.Button('Download RSID Solution Matrix',
                                 className='container',
                                      style={'color':'blue',
                                             'textAlign':'center',
                                             'width':'25%'}),
                        id='download-rsid-solution-link',
                        download="VarCover_RSID_solution_matrix.tsv",
                        href="",
                        target="_blank"),
                html.Hr()
                ])

    return html.Div(div)



# parse uploaded file into table
def parse_vcf(contents, filename, date,
             cost_metric, reduce_singletons):

    if reduce_singletons == 'True':
        reduce_singletons = True
    else:
        reduce_singletons = False

    content_type, content_string = contents.split(',')

    if '.gz' in filename:

        decoded = gzip.decompress(base64.b64decode(content_string))
    else:
        decoded = base64.b64decode(content_string)

    try:
        if 'vcf' in filename:

            user_vcf_path = './user{}.vcf'.format(str(random.random()).split('.')[-1])
            try:
                with open(user_vcf_path, "w+") as f:
                    f.write(decoded.decode("utf-8"))
            except Exception as e:
                print(str(e))

            v = VCF(user_vcf_path,
                    chunksize=1000000,
                    cols=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'FORMAT'])
            v.get_vcf_df_chunk()
            v.add_variant_annotations(inplace=True, drop_hom_ref=False)
            v.df.loc[:, 'GT2'] = v.df['GT2'].fillna('0')
            os.system(' '.join(['rm', user_vcf_path]))

    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
              ])

    df = expand_multiallele(v.df.copy())
    df = create_setcover_df(df)
    vc = varcover(df)
    dropped_vars = vc.dropped_vars.reset_index().rename(columns={0:'GT'})
    soln = vc.getCoverSet(cost=cost_metric,
                          reduceSingletons=reduce_singletons)

    df = reset_index(vc.solution)

    solution_samples = vc.sample_target_allele_cnt.reset_index()

    div = [
            html.H3('VarCover Results for: {}'.format(filename),
                    style={'textAlign': 'left',
                           'color': '#00aeef'}),
            html.H6('Weighting Metric: {}'.format(cost_metric),
                    style={'textAlign': 'left',
                           'color': 'black'}),
            html.H6('Reduce Singletons: {}'.format(reduce_singletons),
                    style={'textAlign': 'left',
                           'color': 'black'}),
            html.H6(datetime.datetime.fromtimestamp(date)),
            html.H4('VarCover Solution Samples',
                    style={'textAlign': 'left',
                           'color': '#d80b8c',
                           'marginBottom':'0.0em'}),
            html.H4('n={}'.format(solution_samples.shape[0]),
                    style={'textAlign': 'left',
                           'color': '#d80b8c',
                           'marginTop':'0.0em'}),
            dt.DataTable(rows=solution_samples.to_dict('records'),
                         filterable=True),
            html.H4('VarCover Solution Matrix',
                    style={'textAlign': 'left',
                       'color': '#d80b8c'}),
            dt.DataTable(rows=df.to_dict('records'), filterable=True)
            ]


    if len(dropped_vars) == 0:

        div.extend([html.H4('No Uncovered Variants',
                            style={'textAlign': 'left',
                                   'color': '#d80b8c'}),
                     html.Br()
                   ])

        return html.Div(div)

    else:
        div.extend([html.H3('Uncovered Variants',
                            style={'textAlign': 'left',
                                   'color': '#d80b8c'}),
                    dt.DataTable(rows=dropped_vars.to_dict('records'),
                                 filterable=True),
                    html.Br()

                   ])


    div.extend([html.A(
                    html.Button('Download VCF Solution Samples',
                                className='container',
                                style={'color':'blue',
                                       'textAlign':'center',
                                       'width':'25%'}),
                        id='download-vcf-sample-link',
                        download="VarCover_VCF_sample_set.tsv",
                        href="",
                        target="_blank",
                        style={'font-size':24,
                              'textAlign':'right'}),

                html.Br(),

                html.A(
                     html.Button('Download VCF Solution Matrix',
                                      className='container',
                                      style={'color':'blue',
                                             'textAlign':'center',
                                             'width':'25%'}),
                        id='download-vcf-link',
                        download="VarCover_VCF_solution_matrix.tsv",
                        href="",
                        target="_blank"),
                html.Hr()
                ])

    return html.Div(div)
