import os,sys,time
import argparse
from varcover import *
from varcover_preprocess import *
from pandasvcf import *



parser = argparse.ArgumentParser()
parser.add_argument('-v','--version', action='version', version='1.0.0')

parser.add_argument('vcf', help='path to vcf file')

parser.add_argument('-c', '--chunksize', default=50000, type=int,
                    help='specify pandasvcf chunksize')

# parser.add_argument('-b','--bed', default='', type=str,
#                    help='restrict to regions specified by bed file')

parser.add_argument('-s','--singletonReduction', action='store_true',
                    help='use singleton reduction modification')

parser.add_argument('-w','--weight', default='standard', choices=['standard','logit'],
                    help='use logit weighting')

parser.add_argument('-n','--niters', default=20, type=int,
                    help='number of PySetCover interations')

parser.add_argument('-o','--output_dir', default='.', type=str,
                    help='specifies the output directory')



if __name__ == '__main__':
    args = parser.parse_args()


    v = VCF(args.vcf, cols=['#CHROM', 'POS', 'REF', 'ID', 'ALT', 'QUAL', 'INFO','FORMAT'], \
         chunksize=args.chunksize)

    v.get_vcf_df_chunk()

    v.add_variant_annotations(inplace=True, drop_hom_ref=False)

    v_df = expand_multiallele(v.df.copy())

    v_df = create_setcover_df(v_df)

    vcover = varcover(v_df)
    vcover.getCoverSet(cost=args.weight,
                       maxit=args.niters,
                       reduceSingletons=args.singletonReduction)

    solution_df = vcover.solution
    print('{}: Samples required to cover {} alleles'.format(solution_df.shape[1],
                                                            solution_df.shape[0]))

    time_stamp = str(time.time())

    if os.path.isdir(args.output_dir):
        pass
    else:
        os.makedirs(args.output_dir)

    write_df_path = os.path.join(args.output_dir, 'varcover_solution_df_{}.tsv'.format(time_stamp))
    solution_df.to_csv(write_df_path,
                       sep='\t')
    print()

    samples = []
    for s in solution_df.columns.unique():
        print(s)
        samples.append(s)

    write_samples_path = os.path.join(args.output_dir, 'varcover_solution_samples_{}.tsv'.format(time_stamp))
    write_sample_f = open(write_samples_path, 'w')
    write_sample_f.write('\n'.join(samples))
    write_sample_f.close()

    allele_cover = os.path.join(args.output_dir, 'varcover_solution_allele_cover_{}.tsv'.format(time_stamp))

    allele_cover_df = solution_df.sum(axis=1).reset_index()
    allele_cover_df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'ALLELE_COUNT']
    allele_cover_df.to_csv(allele_cover,
                           sep='\t',
                           index=False)
