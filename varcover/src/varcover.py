import sys, os
import numpy as np
import pandas as pd
from SetCoverPy import *


class varcover(object):
    '''varcover yields a min-cover set using variant_id x allele count call matrix

    Creates a covering set as self.solution using SetCoverPy

    ToDo:
       1) Implement calculate cost
       2) Identify fixed alt alleles (zero hom-ref)
       3) Implement multiple runs to generate set with best enrichment of rare alleles
       4) Implement report
           4a) Implement clustermap for solution
           4b) Implement summary stats regarding number of rare alleles
           4c) Number of samples
           4d) Variants not covered

    '''
    def __init__(self, df):

        self.df = df
        self.df = self.dropMissingVars()
        self.getAF()
        self.df = self.dropMissingSamples(self.df)


    def dropMissingVars(self):
        '''This function drops variants with no coverage'''
        missing = self.df.sum(axis=1)
        missing = missing[missing<1]
        self.dropped_vars = missing
        return self.df.drop(missing.index, axis=0)


    def getAF(self):
        '''Calculate allele frequencies'''
        af = self.df.sum(axis=1)/(self.df.shape[1]*2)
        af.name = 'af'
        self.af = af


    def dropMissingSamples(self, df):
        '''This function drops Samples with no variants of interest'''
        missing = df.sum(axis=0)
        missing = missing[missing<1]
        self.dropped_samples= missing
        return df.drop(missing.index, axis=1)


    def setSampleTargetAlleleCount(self):
        '''This function creates a Sample x Target Allele Count df'''

        self.sample_target_allele_cnt = pd.DataFrame(self.solution.T.apply(np.sum, axis=1),
                                                     columns=['Target Allele Count'])
        return

    def setTargetAlleleCount(self):
        '''This function creates a Target Allele Count df'''

        self.target_allele_cnt = pd.DataFrame(self.solution.replace(0,np.NaN) \
                                              .stack().reset_index() \
                                              .groupby(['CHROM', 'POS', 'REF', 'ALT'])['sample_ids'] \
                                              .count().sort_values(ascending=False))
        return


    def getCoverSet(self, cost='standard', maxit=20, reduceSingeltons=True):
        '''Creates a covering set using SetCoverPy

           Sets self.solution to covering set matrix.'''


        if cost == 'standard':
            costs = list(self.df.shape[1]*[1.])
        else:
            costs = [float(i) for i in self._calculateCosts(cost=cost)]
        #print(costs)

        costs = pd.Series(costs)
        costs.index = self.df.columns
        self.costs = costs


        print("Solving Set Cover Now:")


        self.testing = 'testing'

        if reduceSingeltons == True:
            self._reduceBySingletons()
            if len(self.singletons_removed) > 0:
                costs = costs.loc[self.singletons_removed.columns]
                g = setcover.SetCover(self.singletons_removed,
                                      cost=costs.values,
                                      maxiters=maxit)
                solution, time_used = g.SolveSCP()
                self.setcov = g
                g_s_df = pd.Series(g.s)
                g_s_df.index = self.singletons_removed.columns
                g_s_df = pd.DataFrame(g_s_df[g_s_df==True])
            else:  # occurs if reduceSingletons collects all samples in self.solution
                self.solution = self.df.loc[self.solution.index, self.solution.columns]
                self.setTargetAlleleCount()
                self.setSampleTargetAlleleCount()
                return
        else:
            g = setcover.SetCover(self.df, cost=costs.values, maxiters=maxit)

            solution, time_used = g.SolveSCP()
            self.setcov = g
            g_s_df = pd.Series(g.s)
            g_s_df.index = self.df.columns
            g_s_df = pd.DataFrame(g_s_df[g_s_df==True])

        if reduceSingeltons == True and len(self.solution) > 0:
            print(self.solution.head())
            print(self.df.head())
            self.solution = self.solution.join(self.df.loc[:, g_s_df.index],
                                              how='outer').fillna(0)

        else: self.solution = self.df[g_s_df.index]

        self.setTargetAlleleCount()
        self.setSampleTargetAlleleCount()

        #self.solution.columns = self.solution.columns.get_level_values(1)


    def _calculateCosts(self, cost='standard', threshold=1.0):
        '''Calculates an allele frequency informed weighting
           for setcoverpy sample cost

           1 / Sum(allele frequency**power)

        '''
        import scipy

        if cost == 'logit':  #logit
            af = pd.Series(self.df.sum(axis=1)/ \
                          (self.df.shape[1]*4)).apply(lambda x: \
                                                      scipy.special.logit(x))
            af = self.df.multiply(af, axis=0).sum(axis=0)
            af = 1 / (af / af.median())

            return af

        if cost == 'standard':
            pass

        # if type(cost):
        #     if cost == 0:  # mean af for each sample
        #         af = self.df.multiply(self.af,axis=0).sum(axis=0)
        #         return 1 / (af / af.median())
        #
        #     else:
        #         af = (self.df.multiply(self.af \
        #                                      .apply(lambda x: x**(power)),axis=0) \
        #                        .sum(axis=0))
        #         af = 1 / (af / af.median())
        #         return af


    def _reduceBySingletons(self):
        '''Isolates samples with singleton variants and
           constructs setcover input dataframe with remaining variants and samples
           '''
        df = self.df.copy()
        singletons = df[df.apply(sum, axis=1)<2].replace(0, np.NaN).stack()
        singletons_samples = singletons.index.get_level_values('sample_ids')
        singletons_vars = singletons.reset_index('sample_ids').index.unique()
        singleton_and_assoc_vars = df.loc[singletons_vars, singletons_samples]

        singletons_removed = df.drop(singletons_vars) \
                               .drop(singletons_samples, axis=1)

        singletons_removed = self.dropMissingSamples(singletons_removed)
        self.singletons_removed = singletons_removed
        self.solution = singleton_and_assoc_vars

        return
