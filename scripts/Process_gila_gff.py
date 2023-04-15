import sys

out_list = []
with open(sys.argv[1]) as f:
	for line in f:
		gene_id = None
		gene_name = None
		line_split = line.strip().split()
		if line_split[0][0] == "#":
			continue
		if line_split[2] == "transcript":
			chrom = line_split[0]
			start = line_split[3]
			attr = line_split[8].split(";")
			for i in attr:
				k = i.split("=")
				if k[0] == "gene_id":
					gene_id = k[1]
				if k[0] == "Name":
					gene_name = k[1]
			if gene_id is not None and gene_name is not None:
				out_list.append([chrom, start, gene_id, gene_name])

with open(sys.argv[2], "w") as o:
	for record in out_list:
		o.write("{}\t{}\t{}\t{}\n".format(record[0], record[1], record[2], record[3]))