import os
import glob
from pathlib import Path

#make a list of all EMBL files

my_path = os.getcwd()
path_embl = os.path.join(my_path, "*.embl")
files_list = glob.glob(path_embl)

#create a new file to output the GC3 values to

file_object1 = open("GC3_content.csv", "w")
file_object1.write("accession,A3,T3,C3,G3,BP,GC3_content\n")

#for each genome, run through all of the FASTA files and calculate GC3

for fle in files_list:
	accession = Path(fle).stem
	
	filepath = os.path.join("/Users/lovedaylewin/Documents/Bath Year 3/FYP/Python_Scripts/learning_python/filtered_genomes/genome_folders/", accession)
	path_tfa = os.path.join(filepath, "*.tfa")
	genes_list = glob.glob(path_tfa)
	
	third_dict = {"A": 0, "T": 0, "C": 0, "G": 0}
	
	for gene_file in genes_list:

		gene = open(gene_file, "r")
		gene_read = (gene.read())
		gene.close()
		
		sequence = (gene_read.split("\n", 1)[1]).upper()
		
		seq_3 = sequence[2::3]
		for s in seq_3:
			try:
				third_dict[s] = third_dict[s] + 1
			except:
				continue
			
	A3 = third_dict["A"]
	T3 = third_dict["T"]
	C3 = third_dict["C"]
	G3 = third_dict["G"]
	BP = A3 + T3 + C3 + G3
	GC3 = ((C3 + G3)/BP)*100
	GC3 = round(GC3, 3)
	
	file_object1.write(f"{accession},{A3},{T3},{C3},{G3},{BP},{GC3}\n")
	
file_object1.close()


