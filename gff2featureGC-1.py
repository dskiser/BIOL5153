#! /home/dskiser/anaconda3/bin/python3.6

# function to clean up a DNA sequence
def clean_seq(input_seq):
	clean = input_seq.upper()
	clean = clean.replace('N','')
	return clean
	
# function to count a nucleotide
def nuc_count(sequence, base):
	length = len(sequence)
	nuc_count = sequence.count(base)
	return nuc_count
	
# check arguments
import sys
if len(sys.argv) < 3:
	print(sys.argv[0] + ": requires .fasta and .gff files")
	sys.exit()

# read arguments
fsa = sys.argv[1]
gff = sys.argv[2]

# find length of genome
fsa_file = open(fsa, "r")
genome = ''
for line in fsa_file:
	if not line.startswith('>'):
		genome = genome + line.strip()
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
	gff_file = open(gff, "r")
	
	# sum of lengths
	length = 0
	G_content = 0
	C_content = 0
	for line in gff_file:
		data = line.split("	")
		if feature == data[2].strip():
			length+=int(data[5])
	# number of G and C nucleotides
			start_index = int(data[3]) - 1
			strand = clean_seq(genome[start_index:int(data[4])])
			G_content += nuc_count(strand, "G")
			C_content += nuc_count(strand, "C")
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
	
	# print output
	features[0] = "exon"
	template = "{0:15}{1:7} {2:7}{3:15}"
	print(template.format(features[feature_num], feature_length[feature_num], percent_genome[feature_num], GC_content[feature_num]))
	feature_num+=1
	
	# close file
	gff_file.close()



