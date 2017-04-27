#! /home/dskiser/anaconda3/bin/python3.6

# function to clean up a DNA sequence
def clean_seq(input_seq):
	clean = input_seq.upper()
	clean = clean.replace('N','')
	return clean
	
# function to count a nucleotide
def nuc_count(sequence, base_1, base_2):
	length = len(sequence)
	base_1_count = sequence.count(base_1)
	base_2_count = sequence.count(base_2)
	combined_base_count = base_1_count + base_2_count
	raw_percent = (combined_base_count/length) * 100
	percent = round(raw_percent, 2)
	return percent

# function to write complementary strand
def complement(seq):
	complement = (seq.replace("A", "t"))
	complement = (complement.replace("C", "g"))
	complement = (complement.replace("T", "a"))
	complement = (complement.replace("G", "c"))
	return complement.upper()
	
# function to get gene name
def get_exon(gff_line):
	break_1 = line.split('Gene')
	break_2 = break_1[1].split(';')
	exon = break_2[0].strip()
	return exon

# key = feature type, value = concatenation of all sequences of that type
feature_sequences = {}
# key = exon, value = sequence	
exon_sequences = {}
# key = gene, value = concatenation of all exons in that gene
gene_sequences = {}

import sys
import collections

# check arguments
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

gff_file = open(gff, "r")
	
for line in gff_file:
	# gather data on feature_types
	data = line.split("	")
	feature = data[2].strip()
	start_index = int(data[3]) - 1
	strand = clean_seq(genome[start_index:int(data[4])])
	if feature in feature_sequences:
		feature_sequences[feature] = feature_sequences[feature] + strand
	else:
		feature_sequences[feature] = strand
	# create dictionary of exons and their sequences
	data = line.split("\t")
	if data[2].strip() == 'CDS':
		exon = get_exon(line)
		start_index = int(data[3]) - 1
		strand = clean_seq(genome[start_index:int(data[4])])
		if data[6].strip() == '-':
			complement_strand = complement(strand)
			exon_sequences[exon] = complement_strand
		else:
			exon_sequences[exon] = strand	
gff_file.close()

# create dictionary of ordered exons
ordered_exons = collections.OrderedDict(sorted(exon_sequences.items()))

# concatenate sequences of each exon
for exon, sequence in ordered_exons.items():
	parsed_exon = exon.split(' ')
	if '-' in parsed_exon[0]:
		parsed_exon = parsed_exon[0].split('-')
	if parsed_exon[0] in gene_sequences:
		gene_sequences[parsed_exon[0]] += sequence
	else:
		gene_sequences[parsed_exon[0]] = sequence

# output
# print feature type information
for feature_type, sequence in feature_sequences.items():
	GC_percent = nuc_count(sequence, "C", "G")
	feature_percent = (len(sequence)/genome_length) * 100
	formated_feature_percent = "(" + str(round(feature_percent, 1)) + "%)"
	template = "{0:15}{1:7} {2:7}{3:15}"
	print(template.format(feature_type, len(sequence), formated_feature_percent, GC_percent))

# print fasta 
print("\n")
for gene, sequence in gene_sequences.items():
	print(">" + gene)
	print(sequence + "\n")
	




