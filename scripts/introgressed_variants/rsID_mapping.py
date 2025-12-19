import argparse
import gzip
import pickle

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', type=str, required=True, help='Path to input file.')
	parser.add_argument('--dict', type=str, required=True, help='Path to rsID mapping dictionary.')
	parser.add_argument('--chr', type=str, required=True, help='Chromosome for which to complete rsID mapping.')
	parser.add_argument('--output', type=str, required=True, help='Path to output file.')
	return parser.parse_args()

def main():
	args = parse_args()

	with gzip.open(args.dict, 'rb') as f:
		rsID_dict = pickle.load(f)

	output_lines = []
	with open(args.input) as infile:
		for line in infile:
			fields = line.strip().split('\t')
			if fields[0] != args.chr:
				continue
			chr_ = fields[0].replace("chr", "")
			pos = fields[1]
			key = (chr_, pos)
			mapped_rsID = rsID_dict.get(key, '.')
			if len(fields) > 2:
				other_fields = fields[2:]
				output_lines.append(f"{chr_}\t{pos}\t{mapped_rsID}\t" + "\t".join(other_fields))
			else:
				output_lines.append(f"{chr_}\t{pos}\t{mapped_rsID}")

	with open(args.output, 'w') as out:
		out.write("\n".join(output_lines) + "\n")

if __name__ == "__main__":
	main()