import os,sys,time
import subprocess
import numpy as np
import pandas as pd
import requests


import dask.dataframe as dd
from dask.delayed import delayed
from functools import partial
from io import StringIO



class EnsemblVar(object):
    """Convenience class to wrap Ensembl Var API
    """

    def __init__(self, rsids,
                 genome='gr37', genotypes=False):

        self.genome = genome
        self.genotypes = genotypes
        self.rsids = rsids
        self.ensembl_json = self._post_ensemblvar()
        self.ensembl_res = self.extract_1KG_genotypes()


    def _post_ensemblvar(self):
        """Post request using self.rsids, self.genome, and self.genotypes
        """

        if self.genome == 'gr37':  # determines genome build to query
            server = "http://grch37.rest.ensembl.org"
        else:
            server = "http://rest.ensembl.org"

        if self.genotypes:  # determines if genotypes should be retrieved
            ext = "/variation/homo_sapiens?genotypes=1"
        else:
            ext = "/variation/homo_sapiens?genotypes=0"

        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        rsid_str = ''.join(['["', '","'.join(self.rsids), '"]'])
        r = requests.post(server+ext, headers=headers, data='{ "ids" : ' + rsid_str +' }')

        if not r.ok:
          r.raise_for_status()
          sys.exit()

        decoded = r.json()
        return decoded


    def extract_1KG_genotypes(self):
        """Construct genotype dataframe from ensembl variant api json
        """

        ensembl_json = self.ensembl_json
        genos = []
        for rsid in ensembl_json.keys():
            mappings = pd.DataFrame(ensembl_json[rsid]['mappings'])
            mappings.index = [rsid] * mappings.shape[0]
            minor_allele = ensembl_json[rsid]['minor_allele']

            if self.genotypes:
                df_geno = pd.DataFrame(ensembl_json[rsid]['genotypes'])
                df_geno = df_geno[df_geno['sample'].str.contains('1000G')]
                df_geno.loc[:, 'sampleid'] = df_geno['sample'].str.split(':').str[-1]
                df_geno.loc[:, 'minor_allele'] = minor_allele
                df_geno.loc[:, 'rsid'] = rsid

                df_geno = df_geno.join(mappings[['allele_string', 'assembly_name', 'location', 'strand']], on='rsid')

                genos.append(df_geno)
            else:
                mappings = mappings[['allele_string', 'assembly_name', 'location', 'strand']]
                mappings.loc[:, 'rsid'] = rsid
                genos.append(mappings[['allele_string', 'assembly_name', 'location', 'strand', 'rsid']])

        return pd.concat(genos, sort=False)


    def set_rsid_bed(self,
                     write_path=''):
        """Get bed file for rsids
        """

        def _parse_location(loc):
            """Parse ensembl location into list
            """
            loc = loc.strip().split(':')
            list_loc = [loc[0]]
            pos = sorted([int(l) for l in loc[1].split('-')])
            list_loc.extend([str(p) for p in pos])

            return list_loc

        def _chrom_int(c):
            try:
                return int(c)
            except:
                if 'X' in c:
                    return 23
                if 'Y' in c:
                    return 24
                return 25

        rsid_bed = pd.DataFrame(self.ensembl_res.set_index('rsid')['location']
                                                 .apply(lambda x: "\t".join(_parse_location(x)))
                                                 .drop_duplicates())
        rsid_bed.columns = ['rsid_bed']
        rsid_bed.loc[:, 'CHROM'] = rsid_bed['rsid_bed'].str.split('\t').str[0]
        rsid_bed.loc[:, 'CHROM_int'] = rsid_bed['CHROM'].apply(_chrom_int)
        rsid_bed.loc[:, 'START'] = rsid_bed['rsid_bed'].str.split('\t').str[1].astype(int) - 1
        rsid_bed.loc[:, 'END'] = rsid_bed['rsid_bed'].str.split('\t').str[-1].astype(int)
        rsid_bed.loc[:, 'rsid_bed'] = rsid_bed['rsid_bed'] + '\t' +rsid_bed.index
        rsid_bed = rsid_bed.sort_values(['CHROM_int', 'START', 'END']).reset_index()

        if write_path != '':
            rsid_bed.reset_index()[['CHROM', 'START', 'END', 'rsid']].to_csv(write_path,
                                                                             sep='\t',
                                                                             header=False,
                                                                             index=False)

        self.rsid_bed = rsid_bed
        return


class DaskBGT(object):
    """Retrive variants intersecting bed file in bgt_dir files
    """
    def __init__(self, bgt_dir, bed,
                bgt_bin='/home/ubuntu/git/bgt/bgt'):

        self.set_bgts_from_dir(bgt_dir)
        self.bed = bed
        self.bgt_bin = bgt_bin



    def get_setcover_df(self):
        """Returns the variant x sample allele count matrix (varcover input)
        """
        self.set_all_bgt_genotypes()
        return self.create_setcover_df()


    def set_bgts_from_dir(self, bgt_dir):
        """Sets self.bcfs by finding all bgt files in bcf_dir
        """
        bgts = [bgt_dir + b.rstrip('.pbf') for b in os.listdir(bgt_dir) if b.endswith('.pbf')]
        self.bgts = sorted(bgts)
        return


    def set_sample_bgt_genotypes(self, sampleID):

        def getGT(file_path, args=""):

            def _cleanBGT(df0, sID):

                df0 = df0[~df0[sID].isin(["2/2", "0/2", "2/0", "1/2", "2/1"])]

                df0.loc[:,"ALT"] = df0["ALT"].str.rstrip(",<M>")
                df0 = df0[~df0.ALT.str.contains('<CN')]
                return df0

            sampleID = args.split(',')[1].split(' ')[0]
            #return sampleID
            arg = [self.bgt_bin, 'view', args, file_path]
            arg = " ".join(arg)
            p0 = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE)
            data = p0.communicate()[0]
            try:
                df = pd.read_table(StringIO(data.decode('utf-8')),\
                                     sep='\t', comment="#", engine="c", header=None)
                df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER", "QUALITY", "INFO", "FORMAT", sampleID]
                df = _cleanBGT(df, sampleID)
                return df

            except:
                return pd.DataFrame()

        daskGT =  partial(getGT, args="-s,{} -B {}".format(sampleID, self.bed))

        dfs = [delayed(daskGT)(f) for f in self.bgts]
        res = dd.compute(dfs, scheduler='processes')
        res = [d for d in res[0] if len(d)>0]
        res = pd.concat(res, sort=False)
        self.bgt_df = res.sort_values(['CHROM', 'POS'])
        return


    def set_all_bgt_genotypes(self):

        def getGTs(file_path, args=""):  # all samples

            def cleanBGT(df0):

                df0 = df0.set_index(df0.columns[:9].tolist()).stack()
                df0 = df0[~df0.isin(["2/2", "0/2", "2/0", "1/2", "2/1"])]
                df0 = df0.unstack().reset_index()

                df0.loc[:,"ALT"] = df0["ALT"].str.rstrip(",<M>")
                df0 = df0[~df0.ALT.str.contains('<CN')]
                return df0

            def _get_sampleids():

                arg = [self.bgt_bin, 'view', args, file_path, '| grep "CHROM" ']
                arg = " ".join(arg)
                p0 = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE)
                sampleids = p0.communicate()[0].decode('utf-8').rstrip('\n').split('\t')[9:]
                return sampleids

            sampleids = _get_sampleids()
            arg = [self.bgt_bin, 'view', args, file_path]
            arg = " ".join(arg)
            print(arg)
            p0 = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE)
            data = p0.communicate()[0]
            try:
                df = pd.read_table(StringIO(data.decode('utf-8')),\
                                     sep='\t', comment="#", engine="c", header=None)
                #sampleids = list(df.columns[9:])
                df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER", "QUALITY", "INFO", "FORMAT"] + sampleids

                df = cleanBGT(df)
                return df

            except:
                return pd.DataFrame()


        daskGT =  partial(getGTs, args="-B {}".format(self.bed))

        dfs = [delayed(daskGT)(f) for f in sorted(self.bgts)]
        res = dd.compute(dfs, scheduler='processes')
        res = [d for d in res[0] if len(d)>0]
        if len(res) == 0:
            assert False
        res = pd.concat(res, sort=False)
        try:
            self.bgt_df = res.sort_values(['CHROM', 'POS'])
        except KeyError:
            self.res = res
        return


    def create_setcover_df(self):
        """Creates variant x sample allele count matrix
        """
        v_df = self.bgt_df.set_index(list(self.bgt_df.columns[:9]))
        v_df = pd.DataFrame(v_df.stack(), columns=['GT'])
        v_df.loc[:, 'GTsum'] = v_df['GT'].apply(lambda x: str(x).count('1'))
        v_df.index.names = list(v_df.index.names)[:-1] + ['sample_ids']
        if 'ID' in v_df.index.names:
            index_cols = ['CHROM', 'POS', 'REF', 'ALT', 'ID']
        else:
            index_cols = ['CHROM', 'POS', 'REF', 'ALT']
        v_df = v_df.pivot_table(index=index_cols,
                                     columns='sample_ids',
                                     values='GTsum',
                                     fill_value=0)

        return v_df
