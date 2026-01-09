import argparse
import gzip
import os
import pickle

from itertools import groupby
from operator import itemgetter

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', type=str, required=True, help='Path to input GWAS summary statistics file.')
	parser.add_argument('--rsID_mapping_directory', type=str, required=True, help='Path to directory with rsID dictionaries.')
	parser.add_argument('--chr_column_index', type=int, required=True, help='Index of chromosome column.')
	parser.add_argument('--pos_column_index', type=int, required=True, help='Index of position column.')
	parser.add_argument('--output', type=str, required=True, help='Path to output file.')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	with open_file(args.input) as input, gzip.open(args.output, 'wt') as out:
		header = input.readline().strip().split('\t')
		out.write('\t'.join(header + ['rsID']) + '\n')

		lines = (line.strip().split('\t') for line in input)
		sorted_lines = sorted(lines, key=itemgetter(args.chr_column_index))

		mapping_cache = {}

		for chr_, group in groupby(sorted_lines, key=itemgetter(args.chr_column_index)):
			chr_str = str(chr_)
			if chr_str == "23":
				chr_str = "X"
			if chr_str not in mapping_cache:
				mapping_file = os.path.join(args.rsID_mapping_directory, f"gnomad_chr{chr_str}_rsID.pkl.gz")
				if not os.path.exists(mapping_file):
					raise FileNotFoundError(f"rsID mapping file not found: {mapping_file}")
				mapping_cache[chr_str] = load_mapping(mapping_file)
			rsID_dict = mapping_cache[chr_str]

			for fields in group:
				while len(fields) < 9:
					fields.append('.')
				if fields[8] == '':
					fields[8] = '.'

				pos_str = str(fields[args.pos_column_index])
				rsID = rsID_dict.get((chr_str, pos_str), '.')

				if len(fields) < 10:
					fields.append(rsID)
				else:
					fields[9] = rsID

				out.write('\t'.join(fields) + '\n')

def open_file(file_path):
	"""Open a file, handling gzipped or plain text files automatically."""
	return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')

def load_mapping(mapping_file):
	"""Load the mapping file into a dictionary with (chr, pos) as keys and rsID as values."""
	open_func = gzip.open if mapping_file.endswith(".gz") else open
	with open_func(mapping_file, 'rb') as file:
		return pickle.load(file)

if __name__ == '__main__':
	main()