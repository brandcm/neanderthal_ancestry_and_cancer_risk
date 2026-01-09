import argparse
import gzip

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', type=str, help='Input decompressed summary stats file.')
	parser.add_argument('--output', type=str, help='Output summary stats file.')
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	with open(args.input, 'rt') as input, gzip.open(args.output, 'wt') as out:
		for i, line in enumerate(input):
			fields = line.strip().split('\t')
			if i == 0:
				out.write('\t'.join(fields + ['rsID']) + '\n')
				continue

			variant_ID = fields[7]
			rsID = variant_ID.split(':')[0]
			rsID = rsID if rsID.startswith('rs') else 'NA'
			out.write('\t'.join(fields + [rsID]) + '\n')

if __name__ == '__main__':
	main()