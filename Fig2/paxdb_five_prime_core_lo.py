'''
purpose of script: read in the 5' and core sequences for each PaxDB genome
-> calculate log odds ratios for each codon's enrichment at the 5'-ends
'''

import math
import os
from pathlib import Path
import glob

#create list of amino acids and dictionary of synonymous codon blocks

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y"]

codons_dict = {
"A": ["GCT", "GCC", "GCA", "GCG"],
"C": ["TGT", "TGC"],
"D": ["GAT", "GAC"],
"E": ["GAA", "GAG"],
"F": ["TTT", "TTC"],
"G": ["GGT", "GGC", "GGA", "GGG"],
"H": ["CAT", "CAC"],
"I": ["ATA", "ATT", "ATC"],
"K": ["AAA", "AAG"],
"L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
"N": ["AAT", "AAC"],
"P": ["CCT", "CCC", "CCA", "CCG"],
"Q": ["CAA", "CAG"],
"R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
"S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
"T": ["ACT", "ACC", "ACA", "ACG"],
"V": ["GTT", "GTC", "GTA", "GTG"],
"Y": ["TAT", "TAC"]
}

#establish file path

my_path = os.getcwd()

#extract all .embl files and use to retrieve accession numbers

path_gbff = os.path.join(my_path, "*.csv")
files_list = glob.glob(path_gbff)

for fle in files_list:
	accession = Path(fle).stem
	print(accession)

#use accession to extract the two .txt files from each accession folder

	filepath = os.path.join("/Users/lovedaylewin/Documents/Loveday_FYP/pax_db_genomes/", accession)
	path_txt = os.path.join(filepath, "*.txt")
	files = glob.glob(path_txt)
	print(files)

#open new output file in each accession folder for log odds csv

	filename = "log_odds_five_prime_core.csv"
	complete_path = os.path.join(filepath, filename)
	file_object1 = open(complete_path, "w")
	file_object1.write(f"codon,amino_acid,{accession}\n")

#read in five prime and core txt files and convert to codons

	five_prime_path = os.path.join(filepath, "five_prime.txt")
	five_prime = open(five_prime_path, "r")
	five_prime_read = five_prime.read()
	n = 3
	five_prime_codons = [five_prime_read[i:i+n] for i in range(0, len(five_prime_read), n)]
	
	core_path = core_path = os.path.join(filepath, "core.txt")
	core = open(core_path, "r")
	core_read = core.read()
	n = 3
	core_codons = [core_read[i:i+n] for i in range(0, len(core_read), n)]
	
#calculate log odds ratio for the enrichment of each codon at the 5' end compared to synonyms
	
	for aa in amino_acids:
		codons = codons_dict[aa]
		alpha = codons
		for c in alpha:
			focal_five_prime = five_prime_codons.count(c)
			focal_core = core_codons.count(c)
			not_c = [codon for codon in codons if not codon.startswith(c)]
			not_c_five_prime = 0
			not_c_core = 0
			for nc in not_c:
				not_c_five_prime = not_c_five_prime + five_prime_codons.count(nc)
				not_c_core = not_c_core + core_codons.count(nc)

#added in if... else as one of the genomes was very GC rich and had a focal_core count of 0
	
			if focal_five_prime and focal_core and not_c_five_prime and not_c_core != 0:
				odds_ratio = (focal_five_prime / focal_core) / (not_c_five_prime / not_c_core)
				log_odds = round((math.log(odds_ratio)), 3)
			else: 
				log_odds = "NA"
			
			file_object1.write(f"{c},{aa},{log_odds}\n")

	file_object1.close()
