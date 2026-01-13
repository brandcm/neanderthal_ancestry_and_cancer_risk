import argparse
import os
import numpy as np
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', type = str, required = True, help = 'Path to input file.')
	parser.add_argument('--assembly', type=str, required=True, help='Human reference assembly ID (hg19 or hg38).')
	parser.add_argument('--chr_column', type=str, required=True, help='Name of chromosome column.')
	parser.add_argument('--pos_column', type=str, required=True, help='Name of position column.')
	parser.add_argument('--rsID_column', type=str, required=True, help='Name of variant ID column.')
	parser.add_argument('--effect_allele_column', type=str, required=True, help='Name of effect allele column.')
	parser.add_argument('--other_allele_column', type=str, required=True, help='Name of other allele column.')
	parser.add_argument('--beta_column', type=str, required=True, help='Name of beta column.')
	parser.add_argument('--p_value_column', type=str, required=True, help='Name of p-value column.')
	parser.add_argument('--output', type = str, required = True, help = 'Path to output file.')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	df = pd.read_csv(
		args.input,
		sep='\t',
		usecols=[args.chr_column, args.pos_column, args.rsID_column, args.effect_allele_column, args.other_allele_column, args.beta_column, args.p_value_column]
	)

	if args.beta_column.lower() in {"odds_ratio", "or"}:
		df["beta"] = np.log(df[args.beta_column])
	else:
		df["beta"] = df[args.beta_column]

	GWAS_SNVs = df[
		(df[args.effect_allele_column].astype(str).str.len() == 1) &
		(df[args.other_allele_column].astype(str).str.len() == 1)
	]

	significant_GWAS_SNVs = GWAS_SNVs[GWAS_SNVs[args.p_value_column] < 5e-8].copy()

	significant_GWAS_SNVs = significant_GWAS_SNVs.rename(columns={
		args.chr_column: 'chr',
		args.pos_column: 'pos',
		args.rsID_column: 'rsID',
		args.effect_allele_column: 'effect_allele',
		args.other_allele_column: 'other_allele',
		args.p_value_column: 'p_value',
	})[['chr', 'pos', 'rsID', 'effect_allele', 'other_allele', 'beta', 'p_value']]

	significant_GWAS_SNVs.to_csv(args.output, sep='\t', index=False, compression='gzip')

if __name__ == '__main__':
	main()