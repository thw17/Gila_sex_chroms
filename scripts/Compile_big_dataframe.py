import argparse
import pandas as pd


def parse_args():
	parser = argparse.ArgumentParser(
		description="This script collects transcript expression data output by "
		"Stringtie and calculates sex-specific mean FPKM values per scaffold.")

	parser.add_argument(
		'--input_files', nargs='+',
		help='Enter the output from Compile_chromstats_results.py '
		'followed by the output from Compile_stringtie_results.py '
		'followed by the output of Calc_het_vcf.py for DNA, '
		'followed by the output of Calc_het_vcf.py for RNA, each'
		'seperated by a space. ')

	parser.add_argument(
		'--output_html_df',
		help='Enter the name you want to use for the html formatted output '
		'dataframe. For example: "gila_big_dataframe.html". ')

	parser.add_argument(
		'--output_csv_df',
		help='Enter the name you want to use for the csv formatted output '
		'dataframe. For example: "gila_big_dataframe.html". ')

	args = parser.parse_args()

	return args


def shifter(df, col_to_shift, pos_to_move):
	# gifted from stackexchange user jpp: https://stackoverflow.com/questions/52616829/how-to-move-a-column-in-a-pandas-dataframe
	arr = df.columns.values
	idx = df.columns.get_loc(col_to_shift)
	if idx == pos_to_move:
		pass
	elif idx > pos_to_move:
		arr[pos_to_move + 1: idx + 1] = arr[pos_to_move: idx]
	else:
		arr[idx: pos_to_move] = arr[idx + 1: pos_to_move + 1]
	arr[pos_to_move] = col_to_shift
	df.columns = arr
	return df


def main():
	args = parse_args()
	input_1 = args.input_files[0]  # Chromstats file
	input_2 = args.input_files[1]  # FPKM data Want to keep the nargs here because it needs multiple files and only the first 2.
	input_3 = args.input_files[2]  # from calc_het_rate.py for DNA
	input_4 = args.input_files[3]  # from calc_het_rate.py for RNA

	df_1 = pd.read_csv(input_1, sep='\t')
	df_2 = pd.read_csv(input_2, sep='\t')
	df_3 = pd.read_csv(input_3, sep='\t')
	df_4 = pd.read_csv(input_4, sep='\t')

	big_df1 = pd.merge(df_1, df_2, on='chrom')
	big_df2 = pd.merge(big_df1, df_3, on='chrom')
	big_df = pd.merge(big_df2, df_4, on='chrom')

	big_df = big_df.drop([
		'G_10_dna.gila2.mkdup.sorted.bam', 'G_16_dna.gila2.mkdup.sorted.bam',
		'G_30_dna.gila2.mkdup.sorted.bam', 'G_35_dna.gila2.mkdup.sorted.bam',
		'G_KI01_dna.gila2.mkdup.sorted.bam', 'G_L_dna.gila2.mkdup.sorted.bam'],
		axis=1)
	big_df = big_df.rename(
		columns={
			'chrom': 'scaffold', 'male_sum': 'male_sum_exp',
			'female_sum': 'female_sum_exp'})

	# big_df = big_df.pipe(shifter, 'length', 1)  # works like a charm

	big_df.to_html(args.output_html_df)
	big_df.to_csv(args.output_csv_df)


if __name__ == "__main__":
	main()
