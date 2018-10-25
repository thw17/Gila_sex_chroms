import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def parse_args():
	parser = argparse.ArgumentParser(
		description="This script takes an input file generated from "
		"chromstats, and creates a dataframe and plots.")

	parser.add_argument(
		'--input_file',
		required=True, help='A text file from chromstats.')

	parser.add_argument(
		'--male_list', nargs='+',
		help='Enter each male sample seperated by a space. '
		'For example: "sample1_dna.ref1.mkdup.sorted.bam sample2_dna.ref1.mkdup.sorted.bam". '
		'These will be fed into a list used to calculate the stats. '
		'Use "--female_list" to enter the female samples.')

	parser.add_argument(
		'--female_list', nargs = '+', help='Enter each female sample seperated by a space. For example: "sample1_dna.ref1.mkdup.sorted.bam sample2_dna.ref1.mkdup.sorted.bam". These will be fed into a list used to calculate the stats. Use "--male_list" to enter the male samples.')

	parser.add_argument(
		'--output_dataframe',required = True, help='Enter the name you would like to use for the csv formatted dataframe. For example: "my_dataframe.txt".')

	parser.add_argument(
		'--output_html', required = True, help='Enter the name you would like to use for the html formatted dataframe. For example: "my_dataframe.html".')

	parser.add_argument(
		'--output_html_no_nan',required = True, help='Enter the name you would like to use for the html formatted dataframe with "nan" & "inf" values omitted. For example: "my_dataframe_no_nan.html".')

	parser.add_argument(
		'--coverage_cutoff', default = 0.95, type = float, help='Enter the cutoff value for the coverage. Values less than the cutoff go into another dataframe (use --output_html_cutoff_df to name that dataframe). ')

	parser.add_argument(
		'--output_html_cutoff_df', required = True, help='Enter the name you would like to use for the html formatted output dataframe holding values less than the value specified by ""--coverage cutoff". For example: "df_ratio_less_than.html". ')

	parser.add_argument(
		'--plot_title', required = True, help='Enter the name of the scatterplot of Female/Male ratio across scaffolds. For example: "gila_chromstats_f_m_ratio.png"'
	)

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	# get data
	inputfile = args.input_file
	print(inputfile)

	df = pd.read_csv(inputfile, sep = '\t')

	# get sum of each column,
	for i in list(df.columns.values):
		if i == 'chrom':
			continue
		Total = df[i].sum()
		df[i] = df[i].apply(lambda x: x / float(Total))

	male_list = args.male_list
	print(male_list)
	female_list = args.female_list

	df['male_average_cov'] = df[male_list].mean(axis=1)
	df['female_average_cov'] = df[female_list].mean(axis=1)
	df['f_m_ratio_cov'] = df.female_average_cov / df.male_average_cov
	df.to_csv(args.output_dataframe, sep='\t', index=False)
	df_no_nan = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)]  # get rid of nans
	# make a df of values that have f_m_ratio less than .95

	df.to_html(args.output_html)
	df_no_nan.to_html(args.output_html_no_nan)

	# if the f_m_ratio is less than .95 add it to another dataframe
	df_ratio_less_than = df_no_nan[
		df_no_nan.f_m_ratio_cov < float(args.coverage_cutoff)]
	df_ratio_less_than.to_html(args.output_html_cutoff_df)

	# ~~~~~~~~~~ plots below ~~~~~~~~~~~~~~
	df_no_nan.plot(kind='scatter', x='chrom', y='f_m_ratio_cov', color='orange')
	plt.title("Female to Male Ratio Across Scaffolds")
	plt.xlabel("Scaffold Number")
	plt.ylabel("Female to Male Ratio")
	plt.savefig(args.plot_title, format="png")


if __name__ == "__main__":
	main()

# ~~~~~~~~~~~~~~~~~ Just for grant to run this thing~~~~~~~~~~~~~~~~~~~~~~~~~

# python gila_chromstats.py --input_file ~/lab_projects/gila1_chrom_stats_count.txt
# --male_list 'G_10_dna.gila1.mkdup.sorted.bam' 'G_16_dna.gila1.mkdup.sorted.bam' 'G_KI01_dna.gila1.mkdup.sorted.bam'
# --female_list 'G_30_dna.gila1.mkdup.sorted.bam' 'G_35_dna.gila1.mkdup.sorted.bam' 'G_L_dna.gila1.mkdup.sorted.bam'
# --output_dataframe 'gila_chromstats_df.txt' --output_html 'gila_chromstats_df.html'
# --output_html_no_nan 'gila_chromstats_df_no_nan.html' --coverage_cutoff '.95'
# --output_html_cutoff_df 'df_ratio_less_than.html' --plot_title 'gila_chromstats_f_m_ratio.png'
