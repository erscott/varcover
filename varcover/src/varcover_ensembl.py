import os,sys,time
import subprocess
import requests
import numpy as np
import pandas as pd
from pandas.io.json import json_normalize


import dask.dataframe as dd
from dask.delayed import delayed
from functools import partial
from io import StringIO



class EnsemblVar(object):
    """Convenience class to wrap Ensembl Var API

    EnsemblVar provides convenience methods for Ensembl Var API and
    methods that reformat results into vcf variant representation.

    Attributes:
        genome: gr37 or gr38

        genotypes: ensembl var api 1000 Genomes Project genotypes

        rsids: iterable of dbSNP ids

        ensembl_res

        ensemb_synonyms

    """

    def __init__(self, rsids,
                 genome='gr37', genotypes=False):
        print('rsid file received')
        self.genome = genome
        self.genotypes = genotypes
        self.rsids = rsids
        self.ensembl_res = self.get_ensembl_res()
        self.ensembl_synonyms = self.add_query_rsid()
        self.multiallelics_split = False
        self.vcf_coords = False
        assert self.ensembl_res['strand'].unique() == 1 #ensures all variants are + strand


        print('ensembl_res generated')

    def get_ensembl_res(self):
        """Retrieves ensembl variant data in batches of 200 variants

        """

        ensembl_res = []
        for rsid_batch in range(0,len(self.rsids), 100):
            print('submitting {}-{} rsids to ensembl api'.format(rsid_batch, rsid_batch+100))
            self.ensembl_json = self._post_ensemblvar(self.rsids[rsid_batch:rsid_batch+100])
            print('Extracting 1KG genotypes')
            ensembl_res.append(self.extract_1KG_genotypes())
        ensembl_res = pd.concat(ensembl_res)
        ensembl_res = ensembl_res.astype('category')
        return ensembl_res


    def _post_ensemblvar(self, rs_batch):
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
        rsid_str = ''.join(['["', '","'.join(rs_batch), '"]'])
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
            mappings = pd.DataFrame(ensembl_json[rsid]['mappings']).astype('category')
            mappings.index = [rsid] * mappings.shape[0]
            minor_allele = ensembl_json[rsid]['minor_allele']
            var_class = ensembl_json[rsid]['var_class']

            if self.genotypes:
                df_geno = pd.DataFrame(ensembl_json[rsid]['genotypes'])
                df_geno = df_geno[df_geno['sample'].str.contains('1000G')].astype('category')
                df_geno.loc[:, 'sample_ids'] = df_geno['sample'].str.split(':').str[-1].astype('category')
                df_geno.loc[:, 'minor_allele'] = minor_allele
                df_geno.loc[:, 'rsid'] = rsid
                df_geno.loc[:, 'var_class'] = var_class

                df_geno = df_geno.join(mappings[['allele_string', 'assembly_name',
                                                 'location', 'strand']],
                                       on='rsid')

                genos.append(df_geno)
            else:
                mappings = mappings[['allele_string', 'assembly_name',
                                     'location', 'strand']]
                mappings.loc[:, 'rsid'] = rsid
                mappings.loc[:, 'var_class'] = var_class
                genos.append(mappings[['allele_string', 'assembly_name',
                                       'location', 'strand', 'rsid',
                                       'var_class']])

        return pd.concat(genos, sort=False).astype('category')


    def split_multiallelic(self):
        """Splits ensembl multiallelic variants into biallelic records

        """

        def _expand_ensembl_multiallele_string(variant):
            """Handles  multiallelic variants into biallelic records
            """

            def expand_snps(variant):
                snp_alleles = []
                ref_base = variant['ref'].unique()[0]
                for allele in variant['allele_string'].unique()[0].split('/'):  #
                    if allele == ref_base:
                        continue
                    allele_variant = variant.copy()
                    allele_variant['alt'] = allele
                    snp_alleles.append(allele_variant)
                return pd.concat(snp_alleles)

            def expand_insertions(variant):
                ins_alleles = []
                for allele in variant['allele_string'].unique()[0] \
                                                      .replace('-/', '') \
                                                      .split('/'):  #
                    allele_variant = variant.copy()
                    allele_variant['alt'] = allele
                    ins_alleles.append(allele_variant)
                return pd.concat(ins_alleles)

            def expand_deletions(variant):
                del_alleles = []
                ref_base = variant['ref'].unique()[0]
                for allele in variant['allele_string'].unique()[0] \
                                                      .replace('-/', '') \
                                                      .split('/'):
                    if allele == ref_base:
                        continue
                    allele_variant = variant.copy()
                    allele_variant['alt'] = allele
                    del_alleles.append(allele_variant)
                return pd.concat(del_alleles)

            split_alleles = []
            var_type = variant['var_class'].unique()[0]
            if var_type == 'SNP':
                split_alleles.append(expand_snps(variant))
            elif var_type == 'insertion':
                split_alleles.append(expand_insertions(variant))
            elif var_type == 'deletion':
                split_alleles.append(expand_deletions(variant))

            return split_alleles

        df = self.ensembl_res.copy()
        self.ensembl_res = self._add_reference_sequence_to_ensembl(df, genome=self.genome)
        del df
        self.ensembl_res.loc[:, 'ref'] = self.ensembl_res['allele_string'].str \
                                              .split('/').str[0]\
                                              .astype('category')  #assumes first allele in allele_string is reference base
        bialleles = self.ensembl_res.groupby('rsid').apply(_expand_ensembl_multiallele_string)
        self.ensembl_res = pd.concat([pd.concat(res) for res in bialleles]).astype('category')
        self.ensembl_res = self.ensembl_res.set_index('query_rsid')
        self.multiallelics_split = True
        return


    def _add_reference_sequence_to_ensembl(self, df, genome='hg19'):
        """Adds a reference sequence column to ensembl_res dataframe
        """

        def get_sequence(chrom, start, end, genome='hg19', just_seq=True):
            '''
            This function requests the nucleotide sequence from the UCSC
            das server(+1 coordinates)
            '''

            def parse_das_response(das):

                seq = [i.decode("utf-8")  for i in das.iter_lines()]
                assert seq[-3] == '</DNA>'
                assert '<DNA length' in seq[4]
                fullseq = ''.join(seq[5:-3])
                assert set(list(fullseq)) - set(['a','t','g','c','N']) == set()
                return fullseq

            if 'chr' not in chrom:
                chrom = 'chr' + str(chrom)

            cmd = ['http://genome.ucsc.edu/cgi-bin/das', genome, \
                   'dna?segment='+chrom+':'+str(start)+','+str(end)]

            das_response = requests.get('/'.join(cmd))

            seq = parse_das_response(das_response)

            if just_seq:
                return seq
            return [i for i in das_response.iter_lines()], seq


        def get_sequence_from_location_prepend(location, genome='hg19'):

            chrom, coords = location.split(':')
            start, end = coords.split('-')

            if int(end) < int(start): # insertions
                prepend_loc = end
                ref = get_sequence(chrom, prepend_loc, prepend_loc,
                                    genome=genome, just_seq=True).upper()
                return ref
            if int(end) > int(start): # insertions
                prepend_loc = int(start) - 1
                ref = get_sequence(chrom, prepend_loc, end,
                                   genome=genome, just_seq=True).upper()
                return ref
            return get_sequence(chrom, start, end,
                                genome=genome, just_seq=True).upper()


        refseq = pd.DataFrame(df['location'].drop_duplicates()) \
                               .set_index('location', drop=False)
        refseq = refseq['location'].apply(get_sequence_from_location_prepend)
        refseq.name = 'vcfref'

        df = df.join(refseq, on='location')
        return df


    def add_vcf_coords(self):
        """Adds vcfCHROM, vcfPOS, vcfREF, vcfalt to dataframe,
           assumes multiallelics_split == True
        """

        if 'alt' in self.ensembl_res and self.multiallelics_split:
            df = self.ensembl_res
        elif self.multiallelics_split == False:
            self.split_multiallelic()
            df = self.ensembl_res

        df_snp = df[df['var_class']=='SNP'].astype('category')
        df_snp.loc[:, 'vcfalt'] = df_snp['alt']
        df_snp.loc[:, 'vcfPOS'] = df_snp['location'].str.split(':').str[1] \
                                      .str.split('-').str[0].astype(int)
        df_snp.loc[:, 'vcfCHROM'] = df_snp['location'].str.split(':').str[0]

        df_ins = df[df['var_class']=='insertion'].astype('category')
        df_ins.loc[:, 'vcfalt'] = df_ins['vcfref'].astype(str) + df_ins['alt'].astype(str)  #adds prepended bases
        df_ins.loc[:, 'vcfPOS'] = df_ins['location'].str.split(':').str[1] \
                                      .str.split('-').str[1].astype(int)
        df_ins.loc[:, 'vcfCHROM'] = df_ins['location'].str.split(':').str[0]

        df_del = df[df['var_class']=='deletion'].astype('category')
        df_del.loc[:, 'vcfalt'] = df_del['vcfref'].str[0]
        df_del.loc[:, 'vcfPOS'] = df_del['location'].str.split(':').str[1] \
                                      .str.split('-').str[0].astype(int) - 1
        df_del.loc[:, 'vcfCHROM'] = df_del['location'].str.split(':').str[0]


        self.ensembl_res = pd.concat([df_snp, df_ins, df_del]).astype('category')
        self.vcf_coords = True
        return

    def add_query_rsid(self):

        rsid_remap = dict()
        for rsid in self.ensembl_json.keys():
            temp_df = self.ensembl_json[rsid]
            synonym_match = set(temp_df['synonyms']) & set(self.rsids)
            if len(synonym_match) > 0:
                rsid_remap[rsid] = list(synonym_match)[0]
        rsid_remap = pd.Series(rsid_remap)
        rsid_remap.name = 'query_rsid'
        self.rsid_remap = rsid_remap
        self.ensembl_res = self.ensembl_res.set_index('rsid') \
                                           .join(rsid_remap, how='left') \
                                           .reset_index() \
                                           .rename(columns={'index':'rsid'})
        self.ensembl_res.loc[self.ensembl_res['query_rsid'].isna(), 'query_rsid'] = \
            self.ensembl_res.loc[self.ensembl_res['query_rsid'].isna(), 'rsid'].values

        return pd.DataFrame(rsid_remap)


    def get_set_cover_df(self):
        """Returns a variant by sample dataframe
        """

        def _get_GTsum_alt(alt_genotype):
            """Counts alt alleles in genotype string

            """
            alt, genotype = alt_genotype
            return genotype.count(alt)

        if not self.genotypes:
            print('No genotypes detected, setting genotypes argument = True')
            self.genotypes = True
            self.ensembl_res = self.get_ensembl_res()
            self.ensembl_synonyms = self.add_query_rsid()
            self.multiallelics_split = False
            self.vcf_coords = False
            assert self.ensembl_res['strand'].unique() == 1 #ensures all variants are + strand

        if not self.multiallelics_split:
            self.split_multiallelic()

        if not self.vcf_coords:
            self.add_vcf_coords()

        self.ensembl_res.loc[:, 'GTsum'] = [_get_GTsum_alt(i) for i in \
                                            self.ensembl_res[['alt', 'genotype']].values]

        self.ensembl_res.rename(columns={'vcfCHROM':'CHROM',
                                         'vcfPOS':'POS',
                                         'vcfref':'REF',
                                         'vcfalt':'ALT'},
                                inplace=True)

        return self.ensembl_res.groupby(['rsid', 'query_rsid','var_class',\
                                        'CHROM', 'POS', 'REF', 'ALT',\
                                         'sample_ids'])['GTsum'].sum().unstack().fillna(0)


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
            rsid_bed.reset_index()[['CHROM', 'START', 'END',
                                    'rsid']] \
                    .to_csv(write_path,
                            sep='\t',
                            header=False,
                            index=False)

        self.rsid_bed = rsid_bed
        return
