#! /home/dskiser/anaconda3/bin/python3.6

# find length of genome
fsa_file = open("watermelon.fsa", "r")
for line in fsa_file:
	if not line.startswith('>'):
		genome = line
fsa_file.close()
genome_length = len(genome)

# create lists of features, feature lengths, percents of genomes, GC
# contents
features = ["CDS","intron","misc_feature","rRNA","repeat_region", "tRNA"]
feature_length = []
percent_genome = []
GC_number = []
GC_content = []

# for each feature, calculate ...
feature_num = 0
for feature in features:
	gff_file = open("watermelon.gff", "r")
	
	# sum of lengths
	length = 0
	G_content = 0
	C_content = 0
	for line in gff_file:
		data = line.split("	")
		if feature == data[2].strip():
			length+=int(data[5])
	# number of G and C nucleotides
			end_index = int(data[4]) + 1
			G_content += genome[int(data[3]):end_index].count("G")
			C_content += genome[int(data[3]):end_index].count("C")
			G_and_C = G_content + C_content
	feature_length.append(length)
	GC_number.append(G_and_C)
	
	# percentage of G and C nucleotides
	raw_percent_GC = (GC_number[feature_num]/feature_length[feature_num]) * 100
	rounded_percent_GC = round(raw_percent_GC, 2)
	GC_content.append(rounded_percent_GC)
			
	# length percentage of genome
	raw_percentage = (feature_length[feature_num]/genome_length) * 100
	rounded_percentage = "(" + str(round(raw_percentage, 1)) + "%)"  
	percent_genome.append(rounded_percentage)
	feature_num+=1
	
	# close file
	gff_file.close()
	
# output
features[0] = "exon"
output_lines = (len(features))
template = "{0:15}{1:7} {2:7}{3:15}"
for number in range(0, output_lines):
	print(template.format(features[number], feature_length[number], percent_genome[number], GC_content[number]))


