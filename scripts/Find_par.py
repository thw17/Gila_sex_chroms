import argparse
import collections
import numpy as np

def parse_args():
	parser = argparse.ArgumentParser(
		description="")

	parser.add_argument(
		"--male_list", nargs="+",
		help="Enter each male sample file seperated by a space.")

	parser.add_argument(
		"--female_list", nargs = "+",
		help="Enter each female sample file seperated by a space")

	parser.add_argument(
		"--output", required = True,
		help="Path to and name of desired output file")

	parser.add_argument(
		"--scaffold", required = True,
		help="Scaffold to analyze")

	args = parser.parse_args()

	return args


def main():

	args = parse_args()

	d = collections.OrderedDict()
	for filename in args.male_list:
		print(filename)
		with open(filename, "r") as i:
			for line in i:
				line1 = line.split()
				if line1[0] == args.scaffold:
					if line1[1] in d:
						d[line1[1]][0].append(float(line1[3]))
					else:
						d[line1[1]] = [[float(line1[3])],[]]
	for filename in args.female_list:
		with open(filename, "r") as i:
			for line in i:
				line1 = line.split()
				if line1[0] == args.scaffold:
					if line1[1] in d:
						d[line1[1]][1].append(float(line1[3]))
					else:
						d[line1[1]] = [[],[float(line1[3])]]

	with open(args.output, "w") as w:
		w.write("start\tmale_depth\tfemale_depth\tratio\tlog2_ratio\n")
		for k in d:
			start = k
			male_depth = np.mean(d[k][0])
			female_depth = np.mean(d[k][1])
			if male_depth > 0:
				ratio = female_depth / male_depth
			else:
				print(d[k])
				continue
			log2_ratio = np.log2(ratio)
			w.write(
				"{}\t{}\t{}\t{}\t{}\n".format(
					start, male_depth, female_depth, ratio, log2_ratio))

if __name__ == "__main__":
	main()
