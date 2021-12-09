import sys

ortho = sys.argv[1]
in_result = sys.argv[2]
outfile = sys.argv[3]

coords = {}
with open(ortho, "r") as f:
	for line in f:
		stripped = line.rstrip()
		split = stripped.split()
		if split[0] != "GFF1_scaff":
			chrom = split[2]
			start = split[3]
			coords[(chrom, start)] = (chrom, start)


with open(outfile, "w") as o:

	with open(in_result, "r") as j:
		for line in j:
			stripped = line.rstrip()
			split = stripped.split()
			if stripped[1] == "chr28":
				chrom = stripped[1]
				start = stripped[3]
				c = (chrom, start)
				if c in coords:
					print(c)
					o.write(line)
				else:
					continue
			else:
				o.write(line)
