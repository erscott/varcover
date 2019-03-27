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
import dash_table_experiments as dt

import pandas as pd
import plotly.graph_objs as go


home = str(Path.home())
sys.path.append(home + '/docker/varcover/varcover/src/')
from varcover_getvars import *
from varcover_preprocess import *
from varcover import *

sys.path.append(home + '/git/pandasVCF')
from pandasvcf import *


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.config['suppress_callback_exceptions'] = True

colors = {
          'background': '#111111',
          'text': '#7FDBFF'
          }


app.layout = html.Div(children=[

  ### HEADER 
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
                  'padding-right':'5px'}),
    

    
  ### WELCOME & CAVEATS
    html.H1(children='Hello Welcome to the VarCover Web App!', 
            style={'textAlign': 'center', 
                   'color': '#d80b8c',
                  'marginBottom': '1.5em'}),
    
    html.H5('Important Note: This is a public demonstration site and provides No Privacy Or Security Features.',
            style={'textAlign':'center', 'color':'#00aeef',
                  'marginBottom': '0.0em'}),
    
    html.H6('Consider using the VarCover package for tasks that require greater data privacy or security.',
            style={'textAlign':'center', 'color':'#00aeef',
                  'marginTop': '0.0em',
                  'marginBottom': '5.0em'}),

    
    
  ### INTERACTIVE COMPONENTS
    html.H2('Please select:',
            style={'textAlign':'center', 'color':'black',
                  'padding-left':'0px'}),

  ### DROP DOWN MENU COST METRIC    
    html.Div(children = [html.H3(' Weighting Metric', style={'textAlign': 'center',
                                   'color':'#d80b8c', 'marginBottom': '0.0em'}),  
                        html.Abbr("hover for help", title="Allele frequency logit weighting increases the likelihood of additional alleles in the solution set with ~ similar sample sizes.")],
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
           'marginBottom': '2.0em'}),
    

    
  ### DROP DOWN MENU REDUCE SINGLETONS
    
    html.Div(children = [html.H3('Reduce Singletons Speedup', style={'textAlign': 'center',
                                   'color':'#d80b8c', 'marginBottom': '0.0em'}),  
                        html.Abbr("hover for help", title="Reduces computational complexity by selecting all samples with singleton alleles prior to solving the min-set cover problem. Most helpful for target sets with >200 variants.")],
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
            html.Abbr("file format ex.", title="#1 rsid per line: .tsv, .csv, .txt file\nrs6500\nrs2032582\nrs56116432")],
             
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
              'margin': '0px'
            },
         # Allow multiple files to be uploaded
         multiple=True
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
    
    html.H6('VCF files: Please do not submit files with private health information or other privacy restrictions', 
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
            multiple=True,
            max_size=1000000
        ),
    
    ### Displaying the uploaded results
    
    html.Div(id='output-vcf-upload'),
    html.Div(dt.DataTable(rows=[{}]), style={'display': 'none'}),
    
    
  ### ATTRIBUTIONS AND FOOTER

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
    
    ]                    
    )
  

    
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
            rsids = pd.read_table(
                               io.StringIO(decoded.decode('utf-8')),
                               sep='\t',
                               header=None,
                               comment='#')
            df.columns = ['rsid']
        elif '.txt' in filename:
            # Assume that the user uploaded a CSV file
            rsids = pd.read_table(
                               io.StringIO(decoded.decode('utf-8')),
                               sep='\t',
                               header=None,
                               comment='#')
            df.columns = ['rsid']
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            rsids = pd.read_excel(io.BytesIO(decoded))
            rsids.columns = ['rsid']
            
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    ensembl = EnsemblVar(rsids['rsid'].tolist())
    
    
    rsid_bed = './rsid{}.bed'.format(str(random.random()).split('.')[-1])
    
    ensembl.set_rsid_bed(rsid_bed)
    
    gts = DaskBGT('~/Documents/rawdata/1000g/phase3/hg19/bgt/',
                  rsid_bed)
    gts = gts.get_setcover_df()
    vc = varcover(gts)
    
    vc.getCoverSet(cost=cost_metric,
                   reduceSingletons=reduce_singletons)
    vc_soln = vc.solution.reset_index()
    vc_soln.loc[:, 'CHROM'] = vc_soln['CHROM'].astype(str)
    vc_soln.loc[:, 'POS'] = vc_soln['POS'].astype(int)
    vc_soln = vc_soln.set_index(['CHROM', 'POS'])
    
    ensembl.rsid_bed = ensembl.rsid_bed.rename(columns={'END':'POS'})
    ensembl.rsid_bed.loc[:, 'POS'] = ensembl.rsid_bed['POS'].astype(int)
    ensembl.rsid_bed.loc[:, 'CHROM'] = ensembl.rsid_bed['CHROM'].astype(str)
    ensembl.rsid_bed = ensembl.rsid_bed.set_index(['CHROM', 'POS'])
    
    vc_soln = vc_soln.join(ensembl.rsid_bed[['rsid']], how='left').reset_index()
    
    os.system(' '.join(['rm', rsid_bed]))
    df = vc_soln.sort_values(['CHROM', 'POS'])
    
    missing_rsids = list(set(rsids['rsid'].tolist()) - set(vc_soln.rsid.values))
    
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
            dt.DataTable(rows=solution_samples.to_dict('records'), filterable=True),
        
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
                    html.Button('Download RSID Solution Samples', className='container',
                                style={'color':'blue',
                                       'textAlign':'center',
                                       'width':'25%'}),
                        id='download-rsid-sample-link',
                        download="VarCover_RSID_sample_set.tsv",
                        href="",
                        target="_blank"),

                html.Br(),

                html.A(
                     html.Button('Download RSID Solution Matrix', className='container',
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
            dt.DataTable(rows=solution_samples.to_dict('records'), filterable=True),
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
                    html.Button('Download VCF Solution Samples', className='container',
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
                     html.Button('Download VCF Solution Matrix', className='container',
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
        children = [
            parse_contents(c, n, d, cost_metric, reduce_singletons) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children
    
    else:
        
        children = html.Div([
            html.H5('Awaiting RSID Upload', 
                    style={'textAlign': 'center', 
                           'color': 'blue'}),
            html.H6('(Max File Size is 1MB)', 
                    style={'textAlign': 'center', 
                           'color': 'blue'})])
    
        return children
    

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
        children = [
            parse_vcf(c, n, d, cost_metric, reduce_singletons) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        
        return children
    
    else:
        
        children = html.Div([
            html.H5('Awaiting VCF Upload', 
                    style={'textAlign': 'center', 
                           'color': 'blue'}),
            html.H6('(Max File Size is 1MB)', 
                    style={'textAlign': 'center', 
                           'color': 'blue'})])
    
        return children

    
@app.callback(
    Output('download-rsid-solution-link', 'href'),
    [Input('output-data-upload', 'children')])
def update_download_link(results):
    
    def _extract_df(parse_child):
        """Extracts datatable data into pandas df from
        output-data-upload
        """
        data = parse_child[0]['props']['children'][8]['props']['rows']
        df = pd.DataFrame(data)

        if 'rsid' in df.columns:
            std_cols = ['rsid', 'CHROM', 'POS', 'ID', 'REF', 'ALT']
            sampleids = list(set(df.columns) - set(std_cols))
            df = df[std_cols + sampleids]
        else:
            std_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
            sampleids = list(set(df.columns) - set(std_cols))
            df = df[std_cols + sampleids]
        return df
    
    
    if type(results) == list:
        
        dff = _extract_df(results)
        csv_string = dff.to_csv(index=False, encoding='utf-8',sep='\t')
        csv_string = "data:text/tsv;charset=utf-8," + urllib.parse.quote(csv_string)
        return csv_string
   

@app.callback(
    Output('download-rsid-sample-link', 'href'),
    [Input('output-data-upload', 'children')])
def update_download_link(results):
    
    def _extract_df(parse_child):
        """Extracts datatable data into pandas df from
        output-data-upload
        """
        data = parse_child[0]['props']['children'][6]['props']['rows']
        df = pd.DataFrame(data)
        df = df[['sample_ids', 'Target Allele Count']]
        return df

    if type(results) == list:
        
        dff = _extract_df(results)
        csv_string = dff.to_csv(index=False, encoding='utf-8',sep='\t')
        csv_string = "data:text/tsv;charset=utf-8," + urllib.parse.quote(csv_string)
        return csv_string
    
    
@app.callback(
    Output('download-vcf-link', 'href'),
    [Input('output-vcf-upload', 'children')])
def update_download_link(results):
    
    def _extract_df(parse_child):
        """Extracts datatable data into pandas df from
        output-data-upload
        """
        data = parse_child[0]['props']['children'][8]['props']['rows']
        df = pd.DataFrame(data)

        if 'rsid' in df.columns:
            std_cols = ['rsid', 'CHROM', 'POS', 'REF', 'ALT']
            sampleids = list(set(df.columns) - set(std_cols))
            df = df[std_cols + sampleids]
        else:
            std_cols = ['CHROM', 'POS', 'REF', 'ALT']
            sampleids = list(set(df.columns) - set(std_cols))
            df = df[std_cols + sampleids]
        return df
    
    
    if type(results) == list:
        
        dff = _extract_df(results)
        csv_string = dff.to_csv(index=False, encoding='utf-8',sep='\t')
        csv_string = "data:text/tsv;charset=utf-8," + urllib.parse.quote(csv_string)
        return csv_string
         

@app.callback(
    Output('download-vcf-sample-link', 'href'),
    [Input('output-vcf-upload', 'children')])
def update_download_link(results):
    
    def _extract_df(parse_child):
        """Extracts datatable data into pandas df from
        output-data-upload
        """
        data = parse_child[0]['props']['children'][6]['props']['rows']
        df = pd.DataFrame(data)
        df = df[['sample_ids', 'Target Allele Count']]
        return df

    if type(results) == list:
        
        dff = _extract_df(results)
        csv_string = dff.to_csv(index=False, encoding='utf-8',sep='\t')
        csv_string = "data:text/tsv;charset=utf-8," + urllib.parse.quote(csv_string)
        return csv_string
        
        
if __name__ == '__main__':
    app.run_server(debug=False)
    
    