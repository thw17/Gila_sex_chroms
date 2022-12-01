import argparse
import collections
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser(
		description="")

	parser.add_argument(
		"--gff1", required=True,
		help="")

	parser.add_argument(
		"--gff2", required=True,
		help="")

	parser.add_argument(
		"--chroms", nargs="+",
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
			if args.chroms is not None:
				if str(parsed[0]) in args.chroms:
					if parsed[2] == "gene":
						scaff = parsed[0]
						start = parsed[3]
						gene = parsed[8].split(';')[1].split('=')[1].upper()
						if gene != "None":
							d[gene] = [scaff, start]
			else:
				if parsed[2] == "gene":
					scaff = parsed[0]
					start = parsed[3]
					gene = parsed[8].split(';')[1].split('=')[1].upper()
					if gene != "None":
						d[gene] = [scaff, start]

	d2 = collections.OrderedDict()
	with open(args.gff2, "r") as f:
		for line in f:
			if line[0] == "#":
				continue
			parsed = line.strip().split()
			if parsed[2] == "region":
				continue
			scaff = parsed[0]
			start = parsed[3]
			gene = parsed[8].split(';')[1].split('=')[1].upper()
			if parsed[2] == "transcript":
				if gene != "None":
					if gene in d:
						identifier = (gene, scaff, start)
						if identifier in d2:
							print(identifier)
						else:
							d2[identifier] = [gene, scaff, start, "transcript", d[gene][0], d[gene][1]]
			elif parsed[2] == "exon":
				if gene != "None":
					if gene in d:
						identifier = (gene, scaff, start)
						if identifier in d2:
							print(identifier)
						else:
							d2[identifier] = [gene, scaff, start, "exon", d[gene][0], d[gene][1]]

	print(len(d))
	print(d)

	print(len(d2))
	print(d2)

	df = pd.DataFrame.from_dict(d2, orient='index', columns=["gene", "GFF2_scaff", "GFF2_start", "type", "GFF1_scaff", "GFF1_start"])
	df_sorted = df.sort_values(["GFF1_scaff", "GFF1_start"])

	df_sorted.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
	main()
