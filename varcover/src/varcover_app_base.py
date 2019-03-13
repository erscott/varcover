import os,sys,random
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import dash_table_experiments as dt
import base64
import datetime
import io

sys.path.append('/Users/ers/docker/varcover/varcover/src/')
from varcover_getvars import *
from varcover_preprocess import *
from varcover import *

sys.path.append('/Users/ers/git/pandasVCF')
from pandasvcf import *



external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

colors = {
          'background': '#111111',
          'text': '#7FDBFF'
          }

df = pd.read_csv(
    'https://gist.githubusercontent.com/chriddyp/'
    'c78bf172206ce24f77d6363a2d754b59/raw/'
    'c353e8ef842413cae56ae3920b8fd78468aa4cb2/'
    'usa-agricultural-exports-2011.csv')

df1 = pd.read_csv(
    'https://gist.githubusercontent.com/chriddyp/' +
    '5d1ea79569ed194d432e56108a04d188/raw/' +
    'a9f9e8076b837d541398e999dcbac2b2826a81f8/'+
    'gdp-life-exp-2007.csv')



app.layout = html.Div(children=[
    html.H1(children='Hello Welcome to the VarCover Web-Interface!', 
            style={'textAlign': 'center', 
                   'color': 'purple'}),

    
    html.H5('Important Note: This is a public demo site and provides No Privacy Or Security Features',
            style={'textAlign':'center', 'color':'red'}),
    html.H5('This website hosts 1000 Genomes Project Phase 3 data',
            style={'textAlign':'center', 'color':'black'}), 
    html.H3('Please select:',
            style={'textAlign':'center', 'color':'black'}),
    html.H4('1) Genome Build,',
            style={'textAlign':'center', 'color':'black'}),
    html.H4('2) RSID or VCF file,',
            style={'textAlign':'center', 'color':'black'}),
    html.H4('Finally, Click the Submit Button to run the VarCover Analysis',
            style={'textAlign':'center', 'color':'black'}),
        
### DROP DOWN MENU    
    html.H3('Genome Build', style={'textAlign': 'center',
                                   'color':'purple'}),   
    dcc.RadioItems(
      options=[
             {'label': 'GRCh37', 'value': 'gr37'},
             {'label': 'GRCh38', 'value': 'gr38'},
             ],
      value='gr37',
     style={'textAlign':'center',  'fontsize':40}),
    
    
###  Inform the User to Submit RSID or VCF
  
   html.H3('Submit an RSID file', 
            style={'textAlign': 'center', 
                   'color': 'purple'}),
    
   html.H6('One RSID Per Line', 
            style={'textAlign': 'center', 
                   'color': 'black'}),
    
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
        html.Div(dt.DataTable(rows=[{}]), style={'display': 'none'}),
    
    
###  Inform the User of Option Submit RSID OR VCF
    
   html.H3('OR', 
            style={'textAlign': 'center', 
                   'color': 'purple'}),
    
   html.H3('Submit a VCF file', 
            style={'textAlign': 'center', 
                   'color': 'purple'}),
    
   html.H6('VCF files: Please do not submit files with private health information or other privacy restrictions', 
            style={'textAlign': 'center', 
                   'color': 'black'}),
    

    ### UPLOAD VCF FILE     
        dcc.Upload(
            id='upload-vcf',
            children=html.Div([
                'VCF file: ',
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
    
    html.Div(id='output-vcf-upload'),
    html.Div(dt.DataTable(rows=[{}]), style={'display': 'none'}),
    
    
    

])
  
    
   
### FUNCTIONS FOR UPLOADING AND TABLE RENDERING

# parse uploaded file into table
def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                             io.StringIO(decoded.decode('utf-8')),
                             header=None)
            df.columns = ['rsid']
        elif 'tsv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_table(
                               io.StringIO(decoded.decode('utf-8')),
                               sep='\t',
                               header=None)
            df.columns = ['rsid']
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
            df.columns = ['rsid']
            
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    ensembl = EnsemblVar(df['rsid'].tolist())
    
    rsid_bed = './rsid{}.bed'.format(str(random.random()).split('.')[-1])
    
    ensembl.set_rsid_bed(rsid_bed)
    
    gts = DaskBGT('/Users/ers/Documents/rawdata/1000g/phase3/hg19/bgt/',
                  rsid_bed)
    gts = gts.get_setcover_df()
    vc = varcover(gts)
    vc.getCoverSet()
    vc_soln = vc.solution.reset_index()
    vc_soln.loc[:, 'CHROM'] = vc_soln['CHROM'].astype(str)
    vc_soln.loc[:, 'POS'] = vc_soln['POS'].astype(int)
    vc_soln = vc_soln.set_index(['CHROM', 'POS'])
    
    ensembl.rsid_bed = ensembl.rsid_bed.rename(columns={'END':'POS'})
    ensembl.rsid_bed.loc[:, 'POS'] = ensembl.rsid_bed['POS'].astype(int)
    ensembl.rsid_bed.loc[:, 'CHROM'] = ensembl.rsid_bed['CHROM'].astype(str)
    ensembl.rsid_bed = ensembl.rsid_bed.set_index(['CHROM', 'POS'])
    
    vc_soln = vc_soln.join(ensembl.rsid_bed[['rsid']], how='left').reset_index()
    
    subprocess.run(['rm', rsid_bed])
    df = vc_soln.sort_values(['CHROM', 'POS'])
    
    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        # Use the DataTable prototype component:
        # github.com/plotly/dash-table-experiments
        dt.DataTable(rows=df.to_dict('records')),

        html.Hr(),  # horizontal line
    ])


# parse uploaded file into table
def parse_vcf(contents, filename, date):
    content_type, content_string = contents.split(',')

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
            subprocess.run(['rm', user_vcf_path])
            
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
              ])
    
    df = expand_multiallele(v.df.copy())
    df = create_setcover_df(df)
    vc = varcover(df)
    vc.getCoverSet()
    df = reset_index(vc.solution)

    return html.Div([
        html.H5(filename),
        html.H6(datetime.datetime.fromtimestamp(date)),

        # Use the DataTable prototype component:
        # github.com/plotly/dash-table-experiments
        dt.DataTable(rows=df.to_dict('records')),

        html.Hr(),  # horizontal line

        # For debugging, display the raw contents provided by the web browser
        #html.Div('Raw Content'),
        #html.Pre(contents[0:200] + '...', style={
        #    'whiteSpace': 'pre-wrap',
        #    'wordBreak': 'break-all'
        #})
    ])


# handle uploaded rsid file and pass to table parser
@app.callback(Output('output-data-upload', 'children'),
              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'),
               State('upload-data', 'last_modified')])
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children
    

# handle uploaded vcf file and pass to table parser
@app.callback(Output('output-vcf-upload', 'children'),
              [Input('upload-vcf', 'contents')],
              [State('upload-vcf', 'filename'),
               State('upload-vcf', 'last_modified')])
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = [
            parse_vcf(c, n, d) for c, n, d in
            zip(list_of_contents, list_of_names, list_of_dates)]
        return children
    


if __name__ == '__main__':
    app.run_server(debug=True)
    