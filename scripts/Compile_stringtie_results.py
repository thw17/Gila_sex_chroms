import argparse
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script collects transcript expression data output by Stringtie and calculates"
        "Sex-specific mean FPKM values per scaffold.")

    parser.add_argument(
        '--sex', nargs = '+', help='A space separated list of sexes ("male" or "female"). Order must be the same as --input_files.')

    parser.add_argument(
        '--input_files', nargs='+', help = 'A space-separated list of t_data.ctab input files.')

    parser.add_argument(
        '--fai', required=True,
        help="Name of and full path to fai file.")

    parser.add_argument(
        '--output_file', required=True,
        help="Name of and full path to output tab-delimited file. Will overwrite if already exists.")

    args = parser.parse_args()

    for x in args.sex:
        if x not in ['male','female']:
            raise ValueError('Sex not specified correctly, must be either "male" or "female".')

    if len(args.sex) != len(args.input_files):
        raise ValueError("Input lengths do not match sex lengths")

    return args


def main():
    t_ids = {}
    t_data = {}

    args = parse_args()
    for idx, file in enumerate(args.input_files):
        with open(file, 'r') as readfile:
            h = readfile.readline()
            for line in readfile:
                temp_line = line.strip()
                temp_split = temp_line.split()

                tid = temp_split[0]
                chr = temp_split[1]
                fpkm = temp_split[-1]

                if tid not in t_ids:
                    t_ids[tid] = chr

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
        if len(t_data[ts]['f']) < 2:
            print('female: ' , ts)
        m_mean = np.mean(np.nan_to_num(t_data[ts]['m']))
        f_mean = np.mean(np.nan_to_num(t_data[ts]['f']))
        mean_table[ts] = {'m': m_mean, 'f': f_mean}

    out_dict = {}
    for ts in t_ids:
        scaff = t_ids[ts]
        if scaff in out_dict:
            out_dict[scaff]['count'] += 1
            out_dict[scaff]['m_sum'] += mean_table[ts]['m']
            out_dict[scaff]['f_sum'] += mean_table[ts]['f']
        else:
            out_dict[scaff] = {'count': 0, 'm_sum': 0, 'f_sum': 0}
            out_dict[scaff]['count'] += 1
            out_dict[scaff]['m_sum'] += mean_table[ts]['m']
            out_dict[scaff]['f_sum'] += mean_table[ts]['f']
    print("done")


    len_dict = {}
    with open(args.fai,'r') as readfile:
        for line in readfile:
            line_split = line.split()
            len_dict[line_split[0]]=line_split[1]

    len_dict_df = pd.DataFrame.from_dict(
        len_dict, orient='index').reset_index().rename(
            index=str, columns={'index': 'scaffold', 0:'length'})

    print(len(mean_table))
    print(len(t_ids))
    table_df = pd.DataFrame.from_dict(out_dict, orient='index').reset_index()
    table_df = table_df.rename(index=str, columns={'index': 'scaffold','count': 'trans_count', 'm_sum': 'male_sum', 'f_sum': 'female_sum'})
    table_df['male_mean'] = table_df.male_sum / table_df.trans_count
    table_df['female_mean'] = table_df.female_sum / table_df.trans_count
    table_df['f_m_ratio'] = table_df.female_mean / table_df.male_mean
    print('done 2')
    out_df = pd.merge(len_dict_df, table_df, on = 'scaffold')

    out_df.to_csv(args.output_file, sep='\t', index=False)
    # table_df.to_csv(args.output_file, sep='\t', index=False)
    # print('done 3')



if __name__ == "__main__":
    main()
