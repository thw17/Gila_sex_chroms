from __future__ import print_function
from __future__ import division

import argparse
import numpy as np
import pandas as pd

samples = [
	"SRR2889297_galgal5", "SRR2889296_galgal5", "SRR2889295_galgal5",
	"SRR2889293_galgal5", "SRR2889292_galgal5", "SRR2889291_galgal5"]

def parse_args():
	parser = argparse.ArgumentParser(
		description="")

	parser.add_argument(
		"--input_file", required=True,
		help="")

	parser.add_argument(
		"--output_file", required=True,
		help="Name of and full path to output tab-delimited file. "
		"Will overwrite if already exists.")

	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	df_1 = pd.read_csv(args.input_file, sep='\t')
	chroms_to_use = ["chr1", "chr2", "chr3", "chr4", "chr5"]
	bools = df_1.scaffold.isin(chroms_to_use)
	df_2 = df_1[bools]
	d = {}
	for i in samples:
		d[i] = df_2[i].mean()
	for i in samples:
		df_1[i] = df_1[i].apply(lambda x: x / d[i])
	df_1["corrected_male_mean"] = df_1[["SRR2889297_galgal5", "SRR2889296_galgal5", "SRR2889295_galgal5"]].mean(axis=1)
	df_1["corrected_female_mean"] = df_1[["SRR2889293_galgal5", "SRR2889292_galgal5", "SRR2889291_galgal5"]].mean(axis=1)
	df_1["corrected_f_m_ratio"] = df_1["corrected_female_mean"] / df_1["corrected_male_mean"]

	if "transcript_id" in list(df_1.columns):
		df_sorted = df_1.sort_values(["scaffold", "transcript_id"], ascending = (True, True))
	elif "exon_id" in list(df_1.columns):
		df_sorted = df_1.sort_values(["scaffold", "exon_id"], ascending = (True, True))

	df_sorted.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
	main()
