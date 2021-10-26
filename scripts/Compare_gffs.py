import argparse
import pandas as pd

def parse_args():
	args = parse_args()

	parser.add_argument(
		"--input_gff1", required=True,
		help="")

	parser.add_argument(
		"--input_gff2", required=True,
		help="")

	parser.add_argument(
		"--chroms", required=True, nargs="+",
		help="")

	parser.add_argument(
		"--output_file", required=True,
		help="Name of and full path to output tab-delimited file. "
		"Will overwrite if already exists.")

	args = parser.parse_args()
	return args

def main():

	args = parse_args()

	d = {}
	with open(args.gff1, "r") as f:
		for line in f:
			if line[0] == "#":
				continue
			parsed = line.strip().split()
			if str(parsed[0]) in args.chroms:
				if parsed[2] == gene:
					scaff = parsed[0]
					start = parsed[3]
					gene = parsed[8].split(';')[1].split('=')[1]
					d[gene] = [scaff, start]

	with open(args.gff2, "r") as f:
		for line in f:
			if line[0] == "#":
				continue
			parsed = line.strip().split()
				scaff = parsed[0]
				start = parsed[3]
				gene = parsed[8].split(';')[1].split('=')[1]
				if gene in d:
					d[gene].append(scaff)
					d[gene].append(start)

	df = pd.DataFrame.from_dict(d, orient='index', columns=["GFF1_scaff", "GFF1_start", "GFF2_scaff", "GFF2_start"])
	df_sorted = df.sort_values(["GFF1_scaff", "GFF1_start")

	df_sorted.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
	main()
