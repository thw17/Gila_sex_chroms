# note that this is hard coded to the gff for anocar2 (ensembl 107 build)

import sys

gff_a = sys.argv[1]
gff_g = sys.argv[2]
result_a = sys.argv[3]
result_g = sys.argv[4]
ortho = sys.argv[5]
outfile = sys.argv[6]

gene_list = {}
transcripts = {}
a_coords = {}
a_flipped = {}
with open(gff_a, "r") as f:
	for line in f:
		id = None
		gene = None
		if line[0] == "#":
			continue
		stripped = line.rstrip()
		split = stripped.split()
		info = split[8].rstrip().split(';')
		if split[2] == "gene":
			if info[0][0:3] == "ID=":
				id = info[0].split(':')[1]
				if info[1][0:5] == "Name":
					gene = info[1].split('=')[1]
					gene_list[id] = gene
			continue
		if split[2] == "mRNA":
			for i in info:
				if i[0:3] == "ID=":
					id = i.split(':')[1]
				if i[0:3] == "Par":
					gene = i.split(':')[1]
			if id is not None and gene is not None:
				transcripts[id] = gene
				chrom = split[0]
				start = split[3]
				a_coords[(chrom, start)] = id
				a_flipped[id] = (chrom, start)

gene_to_transcript = {}
for key in transcripts:
	if transcripts[key] in gene_to_transcript:
		gene_to_transcript[transcripts[key]].append(key)
	else:
		gene_to_transcript[transcripts[key]] = [key]

orthologs = {}
g = []
g_ex = []
a_ex = []
with open(ortho, "r") as f:
	for line in f:
		stripped = line.rstrip()
		split = stripped.split("\t")
		if split[0] == "Orthogroup":
			continue
		gila = ''.join(split[1].rstrip().split()).split(',')
		anolis = ''.join(split[2].rstrip().split()).split(',')
		if len(gila) != 1:
			g_ex.append(gila)
			continue
		if gila[0] in g:
			g_ex.append(gila[0])
			continue
		if gila[0] in g_ex:
			continue
		if len(anolis) > 1:
			a_genes = []
			for k in anolis:
				if k == '':
					print(line)
					print(anolis)
				a_genes.append(transcripts[k])
			if len(set(a_genes)) != 1:
				a_ex.append(a_genes)
				continue
		for k in anolis:
			if k not in a_ex:
				if k == gene_to_transcript[transcripts[k]][0]:
					orthologs[gila[0]] = k

g_coords = {}
with open(gff_g, "r") as f:
	for line in g:
		stripped = line.rstrip()
		split = stripped.split()
		if split[2] == "transcript":
			info = split[8].split(';')
			id = info[0].split('=')[1]
			print(id)
			chrom = split[0]
			start = split[3]
			if id in orthologs:
				g_coords[(chrom, start)] = id

a_results = {}
with open(result_a, "r") as f:
	for line in f:
		stripped = line.rstrip()
		split = stripped.split()
		chrom = str(split[1])
		start = str(split[2])
		male = str(split[3])
		female = str(split[4])
		if (chrom, start) in a_coords:
			a_results[(chrom, start)] = (male, female)

print(len(transcripts))
print(len(a_coords))
print(len(g_coords))
print(len(a_flipped))
print(len(a_results))
print(len(orthologs))

with open(outfile, "w") as o:
	with open(result_g, "r") as j:
		for idx, line in enumerate(j):
			if idx == 0:
				o.write(line + "\tanc_male\tanc_female")
			else:
				stripped = line.rstrip()
				split = stripped.split()
				chrom = str(split[1])
				start = str(split[2])
				if (chrom, start) in g_coords:
					mf = a_coords[a_flipped[orthologs[g_coords[(chrom, start)]]]]
					m1 = mf[0]
					f1 = mf[1]
					o.write(line + "\t{}\t{}".format(m1, f1))
				else:
					continue

