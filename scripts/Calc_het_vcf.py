from __future__ import print_function
import argparse
import collections
import gzip
import statistics
import pandas as pd


def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--vcf", required=True,
		help="")

	parser.add_argument(
		"--sexes", required=True,
		help="")

	parser.add_argument(
		"--min_sites", required=True, type=int,
		help="")

	parser.add_argument(
		"--min_ind", required=True, type=int,
		help="")

	parser.add_argument(
		"--output_file", required=True,
		help="Path to and name of output file")

	parser.add_argument(
		"--suffix", default="",
		help="Text to add to end of column names (except chrom). Dafault is nothing")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()
	dict1 = collections.OrderedDict()
	sample_dict = collections.OrderedDict()

	with gzip.open(args.vcf, 'rt', encoding='utf-8') as f:
		for line in f:
			# Ignore header
			if line[:2] == "##":
				continue
			record = line.split()
			if record[0] == "#CHROM":
				for idx, i in enumerate(record[9:]):
					sample_dict[str(idx)] = i
				continue

			chrom = record[0]
			if chrom not in dict1:
				dict1[chrom] = []
				for i in range(len(sample_dict)):
					dict1[chrom].append([0, 0])

			# Looking for biallelic snps
			if len(record[3]) != 1:
				continue
			alts = record[4].split(",")
			alts = [x for x in alts if x != "<NON_REF>"]
			if len(alts) not in [1, 2]:
				continue
			lengths = [len(x) for x in alts]
			if len(lengths) == 0:
				continue
			if all(x == 1 for x in lengths) is False:
				continue

			for idx, i in enumerate(record[9:]):
				format = i.split(":")
				gt = format[0]
				if gt == "./.":
					continue
				if gt == "0/0":
					continue
				dict1[chrom][idx][0] += 1
				if gt[0] != gt[2]:
					dict1[chrom][idx][1] += 1

	sex_dict = collections.OrderedDict()
	with open(args.sexes, "r") as f:
		for line in f:
			line_split = line.strip().split()
			if len(line_split) != 2:
				continue
			sex_dict[line_split[0]] = line_split[1]

	for i in sample_dict:
		sample_dict[i] = sex_dict[sample_dict[i]]

	results_dict = collections.OrderedDict()
	for i in dict1:
		male_vals = []
		fem_vals = []

		for idx, k in enumerate(dict1[i]):
			sex = sample_dict[str(idx)]
			if k[0] >= args.min_sites:
				if sex == "male":
					male_vals.append(float(k[1]) / k[0])
				elif sex == "female":
					fem_vals.append(float(k[1]) / k[0])
				else:
					raise ValueError("Sexes need to be the strings 'male' or 'female'.")

		if len(male_vals) >= args.min_ind and len(fem_vals) >= args.min_ind:
			results_dict[i] = [i, statistics.mean(male_vals), statistics.mean(fem_vals)]
		else:
			results_dict[i] = [i, -1, -1]

		print(i)
		print(results_dict[i])

	df = pd.DataFrame.from_dict(results_dict, orient='index')
	df.columns = [
		"chrom",
		"male_het_rate_{}".format(args.suffix),
		"female_het_rate_{}".format(args.suffix)]
	print(df)
	df = df.sort_values(["chrom"])
	df.to_csv(args.output_file, sep='\t', index=False)


if __name__ == "__main__":
	main()
