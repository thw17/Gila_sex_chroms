import sys

ortho = sys.argv[1]
gff1_result = sys.argv[2]
gff2_result = sys.argv[3]
outfile = sys.argv[4]
region_type = sys.argv[5]

gff1_coords = {}
gff2_gff1_lookup = {}
gff1_mult = []
gff2_mult = []
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
				if (chrom_1, start_1) in gff1_coords:
					del gff1_coords[(chrom_1, start_1)]
					gff1_mult.append((chrom_1, start_1))
					continue
				if (chrom_1, start_1) in gff1_mult:
					continue
				gff1_coords[(chrom_1, start_1)] = (chrom_1, start_1)
				if (chrom_2, start_2) in gff2_gff1_lookup:
					del gff2_gff1_lookup[(chrom_2, start_2)]
					gff2_mult.append((chrom_2, start_2))
					del gff1_coords[(chrom_1, start_1)]
					continue
				if (chrom_2, start_2) in gff2_mult:
					continue
				gff2_gff1_lookup[(chrom_2, start_2)] = (chrom_1, start_1)
			else:
				continue

for k in list(gff2_gff1_lookup):
	if gff2_gff1_lookup[k] in gff1_mult:
		gff2_mult.append(k)
		del gff2_gff1_lookup[k]

# gff2_gff1_lookup_final = {key: val for key, val in gff2_gff1_lookup if key not in gff2_mult}

print("Multiple hits in gff1: {}".format(gff1_mult), len(gff1_mult))
print("")
print("Multiple hits in gff2: {}".format(gff2_mult), len(gff2_mult))

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
