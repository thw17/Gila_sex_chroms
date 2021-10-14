

def grouped(iterator, size):
	yield tuple(next(iterator) for _ in range(size))

scaff_dict = {}
with open("komodo_scaff218_gila2_align.maf", "r") as maf:
	block_count = 0
	line_counter = 1
	for line in maf:
		if line[0] == "#":
			continue
		else:
			if line_counter == 1:
				score = line.strip().split()[1].split("=")[1]
			if line_counter == 4:
				line_counter = 1
				block_count += 1
			elif line_counter == 2:
				parsed = line.strip().split()
				scaff = parsed[1]
				length = parsed[3]
				if scaff in scaff_dict:
					scaff_dict[scaff] += int(length)
				else:
					scaff_dict[scaff] = int(length)
				line_counter += 1
				print(len(scaff_dict))
			else:
				line_counter += 1
print("Maf blocks: {}".format(block_count))
print("Scaffold lengths:")
for i in scaff_dict:
	print("{}: {}".format(i, scaff_dict[i]))
