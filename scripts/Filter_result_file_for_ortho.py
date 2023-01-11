import sys

ortho = sys.argv[1]
in_result = sys.argv[2]
outfile = sys.argv[3]
region_type = sys.argv[4]

coords = {}
with open(ortho, "r") as f:
	for line in f:
		stripped = line.rstrip()
		split = stripped.split()
		if split[0] not in ["gene", "NONE"]:
			if split[3] == region_type:
				chrom = split[1]
				start = split[2]
				coords[(chrom, start)] = (chrom, start)
			else:
				continue

with open(outfile, "w") as o:

	with open(in_result, "r") as j:
		for line in j:
			stripped = line.rstrip()
			split = stripped.split()
			if split[1] == "chr28":
				chrom = split[1]
				start = split[2]
				c = (chrom, start)
				if c in coords:
					o.write(line)
				else:
					print(c)
					continue
			else:
				o.write(line)
