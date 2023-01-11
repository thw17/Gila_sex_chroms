import sys

ortho = sys.argv[1]
gff1_result = sys.argv[2]
gff2_result = sys.argv[3]
outfile = sys.argv[4]
region_type = sys.argv[5]

gff1_coords = {}
gff2_gff1_lookup = {}
with open(ortho, "r") as f:
	for line in f:
		stripped = line.rstrip()
		split = stripped.split()
		if split[0] not in ["gene", "NONE"]:
			if split[3] == region_type:
				chrom_1 = split[4]
				start_1 = split[5]
				chrom_2 = split[1]
				start_2 = split[2]
				gff1_coords[(chrom_1, start_1)] = (chrom_1, start_1)
				if (chrom_2, start_2) in gff2_gff1_lookup:
					print(line)
				gff2_gff1_lookup[(chrom_2, start_2)] = (chrom_1, start_1)
			else:
				continue

gff2_in_gff1 = {}
with open(gff2_result, "r") as g:
	for line in g:
		stripped = line.rstrip()
		split = stripped.split()
		chrom = split[1]
		start = split[2]
		male = split[3]
		female = split[4]
		if (chrom, start) in gff2_gff1_lookup:
			gff2_in_gff1[gff2_gff1_lookup[(chrom, start)]] = (male, female)

print("gff1_coords", len(gff1_coords))
print("gff2_gff1_lookup", len(gff2_gff1_lookup))
print("gff2_in_gff1", len(gff2_in_gff1))

with open(outfile, "w") as o:

	with open(gff1_result, "r") as j:
		for idx, line in enumerate(j):
			if idx == 0:
				o.write(line + "\tanc_male\tanc_female")
			else:
				stripped = line.rstrip()
				split = stripped.split()
				chrom = split[1]
				start = split[2]
				if (chrom, start) in gff1_coords:
					m1 = gff2_in_gff1[(chrom, start)][0]
					f1 = gff2_in_gff1[(chrom, start)][1]
					o.write(line + "\t{}\t{}".format(m1, f1))
				else:
					continue
