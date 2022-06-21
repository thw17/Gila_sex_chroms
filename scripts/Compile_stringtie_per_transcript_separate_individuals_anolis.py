from __future__ import print_function
from __future__ import division

import argparse
import numpy as np
import pandas as pd


def parse_args():
	parser = argparse.ArgumentParser(
		description="This script collects transcript expression data output "
		"by Stringtie and calculates"
		"Sex-specific mean FPKM values per scaffold.")

	parser.add_argument(
		"--sex", nargs="+",
		help="A space separated list of sexes ('male' or 'female'). "
		"Order must be the same as --input_files.")

	parser.add_argument(
		"--input_files", nargs="+",
		help="A space-separated list of t_data.ctab input files.")

	parser.add_argument(
		"--suffix", required=True,
		help="String of text to append to the end of relevant column names. "
		"E.g., '--suffix denovo' would add '_denovo' to the end of count, sum, "
		"mean, and ratio columns.")

	parser.add_argument(
		"--output_file", required=True,
		help="Name of and full path to output tab-delimited file. "
		"Will overwrite if already exists.")

	args = parser.parse_args()

	for x in args.sex:
		if x not in ['male', 'female']:
			raise ValueError(
				'Sex not specified correctly, must be either "male" or "female".')

	if len(args.sex) != len(args.input_files):
		raise ValueError("Input lengths do not match sex lengths")

	return args


def main():
	t_ids = {}
	t_data = {}
	sample_order = [[],[]]

	args = parse_args()
	for idx, file in enumerate(args.input_files):
		# This if/else block is hard coded to the expected path from the Snakefile
		if args.sex[idx] == 'male':
			a = file.split("/")[1]
			print(a)
			b = a.split("_rna_")[0]
			print(b)
			sample_order[0].append(b)
		else:
			a = file.split("/")[1]
			b = a.split("_rna_")[0]
			sample_order[1].append(b)

		with open(file, 'r') as readfile:
			h = readfile.readline()
			for line in readfile:
				temp_line = line.strip()
				temp_split = temp_line.split()

				tid = temp_split[0]
				chr = temp_split[1]
				fpkm = temp_split[-1]
				start = temp_split[3]

				if tid not in t_ids:
					t_ids[tid] = [chr,start]

				if tid not in t_data:
					t_data[tid] = {'m': [], 'f': []}
				if args.sex[idx] == 'male':
					t_data[tid]['m'].append(float(fpkm))
				else:
					t_data[tid]['f'].append(float(fpkm))

	mean_table = {}
	for ts in t_data:
		if len(t_data[ts]['m']) < 2:
			print('male: ', ts)
			continue
		if len(t_data[ts]['f']) < 2:
			print('female: ', ts)
			continue
		m_mean = np.mean(np.nan_to_num(t_data[ts]['m']))
		f_mean = np.mean(np.nan_to_num(t_data[ts]['f']))
		chrom = t_ids[ts][0]
		start = t_ids[ts][1]
		mean_table[ts] = {
			'scaffold': chrom, 'start': start, 'male_mean_expn': m_mean, 'female_mean_expn': f_mean,
			sample_order[0][0]: t_data[ts]['m'][0],
			sample_order[0][1]: t_data[ts]['m'][1],
			sample_order[1][0]: t_data[ts]['f'][0],
			sample_order[1][1]: t_data[ts]['f'][1]}

	table_df = pd.DataFrame.from_dict(mean_table, orient='index').reset_index()
	table_df = table_df.rename(
		index=str,
		columns={
			'index': 'transcript_id'})

	table_df[
		'f_m_ratio_{}'.format(
			args.suffix)] = table_df.female_mean_expn / table_df.male_mean_expn
	table_df = table_df.rename(
		index=str,
		columns={
			'trans_count': 'trans_count_{}'.format(args.suffix),
			'male_mean_expn': 'male_mean_{}'.format(args.suffix),
			'female_mean_expn': 'female_mean_{}'.format(args.suffix)})

	table_df.to_csv(args.output_file, sep='\t', index=False)


if __name__ == "__main__":
	main()
