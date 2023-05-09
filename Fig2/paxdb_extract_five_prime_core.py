'''
Aim of script: output a txt file containing a string of all the 5' codons and 
a txt file containing a string of all the core codons for each of the PaxDB genomes
'''

import os
import glob
from pathlib import Path
from Bio import SeqIO

#make a list of the GenBank files 

my_path = os.getcwd()

path_gbff = os.path.join(my_path, "pax_db_genomes", "*.gbff")
files_list = glob.glob(path_gbff)
	
for fle in files_list:

#use BioPython to extract the CDS and their locus tags -> convert sequence into codons
#-> extract first 11 codons of each gene (excluding start codon) and join in a string
#-> extract middle 11 codons of each gene and join in a string -> output to txt files 

	for gb_record in SeqIO.parse(open(fle, "r"), "genbank"):
		print(f"Name: {gb_record.name}, features: {len(gb_record.features)}")
		
		filepath = os.path.join("/Users/lovedaylewin/Documents/Bath Year 3/FYP/Python_Scripts/learning_python/filtered_genomes/pax_db_genomes/", gb_record.name)

#do the B. subtilis separately as the PaxDB abundance data uses the old locus tag

		if gb_record.name == "NC_000964":
		
			filename1 = "five_prime.txt"
			complete_path1 = os.path.join(filepath, filename1)
			
			filename2 = "core.txt"
			complete_path2 = os.path.join(filepath, filename2)
			
			file_object1 = open(complete_path1, "w")
			file_object2 = open(complete_path2, "w")

			for feature in gb_record.features:
				if feature.type == "CDS":
					locus_tag = feature.qualifiers["old_locus_tag"][0]
					gene_sequence = feature.extract(gb_record.seq)
				
					five_prime_string = ""
					core_string = ""
					
					sequence = str(gene_sequence)
					n = 3
					codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]		
					
					five_prime = codons[1:12]
					five_prime_codon_string = "".join(five_prime)
					five_prime_string = five_prime_string + five_prime_codon_string
					
					K = 11
					beg_indx = (len(codons) // 2) - (K // 2)
					end_indx = (len(codons) // 2) + (K // 2)
					core = codons[beg_indx:end_indx + 1]
					core_codon_string = "".join(core)
					core_string = core_string + core_codon_string
					
					file_object1.write(five_prime_string)
					file_object2.write(core_string)
				
			file_object1.close()
			file_object2.close()
				
		else: 
		
			filename1 = "five_prime.txt"
			complete_path1 = os.path.join(filepath, filename1)
			
			filename2 = "core.txt"
			complete_path2 = os.path.join(filepath, filename2)
			
			file_object1 = open(complete_path1, "w")
			file_object2 = open(complete_path2, "w")
	
			for feature in gb_record.features:
				if feature.type == "CDS":
					locus_tag = feature.qualifiers["locus_tag"][0]
					gene_sequence = feature.extract(gb_record.seq)
				
					five_prime_string = ""
					core_string = ""
					
					sequence = str(gene_sequence)
					n = 3
					codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]		
					
					five_prime = codons[1:12]
					five_prime_codon_string = "".join(five_prime)
					five_prime_string = five_prime_string + five_prime_codon_string
					
					K = 11
					beg_indx = (len(codons) // 2) - (K // 2)
					end_indx = (len(codons) // 2) + (K // 2)
					core = codons[beg_indx:end_indx + 1]
					core_codon_string = "".join(core)
					core_string = core_string + core_codon_string
					
					file_object1.write(five_prime_string)
					file_object2.write(core_string)
					
			file_object1.close()
			file_object2.close()
