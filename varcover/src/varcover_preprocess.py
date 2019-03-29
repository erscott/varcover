import sys, os
import pandas as pd

#Import pdVCFsingle package
homePath = os.path.expanduser('~')

PANDAS_VCF_SRC_DIR = homePath + '/git/pandasVCF/'
sys.path.append( PANDAS_VCF_SRC_DIR )
from pandasvcf import *
pd.options.mode.chained_assignment = None #supressing the chained assignment warnings
pd.set_option('display.precision',30)


def clean_alt(ref_a1_a2_phase):
    '''
    This function cleans the ALT and GT columns of multi-allelic
    variants

    '''
    ref, a1, a2, phase = ref_a1_a2_phase
    newAlt = []
    newGT = []
    for n,allele in enumerate([a1,a2]):
        if allele != ref:
            if allele not in newAlt:
                newAlt.append(allele)
                newGT.append(str(len(newAlt)))
                continue
            else:
                newGT.append(str(len(newAlt)))
                continue
        newGT.append('0')
    return (','.join(newAlt), newGT[0], newGT[1])


def expand_multiallele(df):


    def alt_count(df):
        '''Count max number of alternative alleles'''
        return df['ALT'].str.count(',').max()


    def clean_alt(ref_a1_a2_phase):
        '''
        This function cleans the ALT and GT columns of multi-allelic
        variants

        '''
        ref, a1, a2, phase = ref_a1_a2_phase
        newAlt = []
        newGT = []
        for n,allele in enumerate([a1,a2]):
            if allele != ref:
                if allele not in newAlt:
                    newAlt.append(allele)
                    newGT.append(str(len(newAlt)))
                    continue
                else:
                    newGT.append(str(len(newAlt)))
                    continue
            newGT.append('0')
        return ','.join(newAlt), newGT[0], newGT[1]


    def slow_expand_multiallele(line):
        '''This function splits multi-alt allele variants
        into separate rows.  It is slow due to a for loop
        '''
        newlines = []
        alts = str(line['newALT'])
        for n,allele in enumerate(alts.split(',')):
            newline = line.copy()
            assert int(line['newGT2'])>1
            assert int(line['newGT1'])==1
            if n == 0:

                newline['newALT'] = allele
                newline['newGT2'] = 0
                newlines.append(newline)
            else:
                newline['newALT'] = allele
                newline['newGT1'] = 0
                newline['newGT2'] = 1
                newlines.append(newline)

        return pd.DataFrame(newlines)


    def fast_expand_multiallele(ma_df):
        '''This function splits multi-alt allele variants
        into separate rows.  It is fast due to vector assignments
        '''
        alt2_df = ma_df.copy()

        assert int(ma_df['newGT2'].max())==2
        assert int(ma_df['newGT1'].max())==1
        ma_df.loc[:, 'newALT'] = ma_df['newALT'].str.split(',').str[0]
        ma_df.loc[:, 'newGT2'] = 0


        alt2_df.loc[:, 'newALT'] = alt2_df['newALT'].str.split(',').str[1]
        alt2_df.loc[:, 'newGT1'] = 0
        alt2_df.loc[:, 'newGT2'] = 1

        return ma_df.append(alt2_df).sort_index()



    # Prepare the multi-allelic dataframe
    df.loc[:, 'newALT'] = df.reset_index()['ALT'].values
    df.loc[:, 'newGT1'] = df['GT1'].astype(int)
    df.loc[:, 'newGT2'] = df['GT2'].astype(int)
    ma_df = df.query('multiallele>0')
    if len(ma_df) > 0:
        df = df.drop(ma_df.index).reset_index()

        ma_df.reset_index(inplace=True)
        ma_df.loc[:, ['newALT', 'newGT1', 'newGT2']] = [i for i in map(clean_alt,                                                       ma_df[['REF', 'a1', 'a2', 'phase']].values)]
        ma_df = ma_df[ma_df.newALT.str.contains(',')]

        # Triage through fast vector or slow for loop
        if alt_count(ma_df) == 1:  # if only triallelic variants in multiallele df then can use fast vector function

            expanded_df = fast_expand_multiallele(ma_df)

            # Combine results with non-multi-allelic snps
            df = df.append(expanded_df)

        else:  # if only >triallelic variants in multiallele df then must use slow for-loop function

            biallelic_vars = ma_df[ma_df['ALT'].str.count(',')==1]
            multiallelic_vars = ma_df.drop(biallelic_vars.index)

            res = []
            for i in multiallelic_vars.index:  # iteratate through variants
                l = multiallelic_vars.loc[i]
                res.append(slow_expand_multiallele(l))

            try:
                expanded_df = pd.concat(res)
                expanded_df = expanded_df.append(biallelic_vars)

                # Combine results with non-multi-allelic snps
                df = df.append(expanded_df)
            except ValueError: # occurs if no multiallelic_vars
                pass

    else:
        df = df.reset_index()

    return df.sort_values(['CHROM', 'POS', 'sample_ids','newALT'])


def create_setcover_df(df):

    df.loc[:, 'GTsum'] = df['newGT1'].astype(int) + df['newGT2'].astype(int)
    del df['ALT']
    df.rename(columns={'newALT':'ALT'}, inplace=True)
    df = df.drop_duplicates(['CHROM', 'POS', 'REF', 'ALT', 'sample_ids', 'GTsum'])
    return df.set_index(['CHROM', 'POS', 'REF', 'ALT', 'sample_ids'])['GTsum'] \
                    .unstack().fillna(0)


def create_setcover_df_from_vcf(vcf_df):
    """Creates variant x sample allele count matrix using standard vcf dataframe
    """
    v_df = vcf_df.set_index(list(vcf_df.columns[:9]))
    v_df = pd.DataFrame(v_df.stack(), columns=['GT'])
    v_df.loc[:, 'GTsum'] = v_df['GT'].apply(lambda x: str(x).count('1'))
    v_df.index.names = list(v_df.index.names)[:-1] + ['sample_ids']
    v_df = v_df.pivot_table(index=['CHROM', 'POS', 'REF', 'ALT'],
                                 columns='sample_ids',
                                 values='GTsum',
                                 fill_value=0)
    return v_df


def reset_index(df):
    '''Returns DataFrame with index as columns'''
    index_df = df.index.to_frame(index=False)
    df = df.reset_index(drop=True)
    #  In merge is important the order in which you pass the dataframes
    # if the index contains a Categorical.
    # pd.merge(df, index_df, left_index=True, right_index=True) does not work
    return pd.merge(index_df, df, left_index=True, right_index=True)
