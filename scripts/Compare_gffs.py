import argparse
import pandas as pd

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
