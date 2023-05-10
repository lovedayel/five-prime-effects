import os
import glob
from pathlib import Path

#establish file path

my_path = os.getcwd()

#glob for all files ending in .embl and use the resulting list to extract accession numbers

path_embl = os.path.join(my_path, "*.embl")
files_list = glob.glob(path_embl)

for fle in files_list:
	accession = Path(fle).stem
	
#for each accession folder, extract all files that end in .tfa (the gene files)

	filepath = os.path.join("/Users/lovedaylewin/Documents/Bath Year 3/FYP/Python_Scripts/learning_python/filtered_genomes/genome_folders/", accession)
	path_tfa = os.path.join(filepath, "*.tfa")
	genes_list = glob.glob(path_tfa)

#open new txt files to output the five prime and core sequences for each accession

	prime_filename = "five_prime.txt"
	prime_complete_path = os.path.join(filepath, prime_filename)
	file_object1 = open(prime_complete_path, "w")
	
	core_filename = "core.txt"
	core_complete_path = os.path.join(filepath, core_filename)
	file_object2 = open(core_complete_path, "w")

	for gene_file in genes_list:

		gene = open(gene_file, "r")
		gene_read = (gene.read())
		gene.close()
	
#extract the pure sequence from each gene file and convert to codons
	
		sequence = (gene_read.split("\n", 1)[1]).upper()

		five_prime_string = ""
		core_string = ""
	 
		n = 3
		codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]

#extract the first 11 codons (excluding the start codon), and the middle 11 codons

		five_prime = codons[1:12]
		five_prime_codon_string = "".join(five_prime)
		five_prime_string = five_prime_string + five_prime_codon_string
	
		K = 11
		beg_indx = (len(codons) // 2) - (K // 2)
		end_indx = (len(codons) // 2) + (K // 2)
		core = codons[beg_indx:end_indx + 1]
		core_codon_string = "".join(core)
		core_string = core_string + core_codon_string

#write the sequences to the output files
	
		file_object1.write(five_prime_string)
		file_object2.write(core_string)
	
file_object1.close()
file_object2.close()



