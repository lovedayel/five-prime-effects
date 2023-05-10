from Bio import SeqIO
from Bio.SeqUtils import GC123
import os
import glob
from pathlib import Path

#make a list of the Genbank annotations for the five PaxDb genomes

my_path = os.getcwd()

path_gbff = os.path.join(my_path, "pax_db_genomes", "*.gbff")
files_list = glob.glob(path_gbff)

file_object2 = open("paxdb_gc3.csv", "w")
file_object2.write("accession,A3,T3,C3,G3,BP,GC3_content\n")
	
for fle in files_list:

#for each genome, extract each CDS by its locus tag and extract the first 11 codons
#output a csv for each genome containing the 5' sequences
#calculate GC3 for each bacteria and output results to a second csv

	for gb_record in SeqIO.parse(open(fle, "r"), "genbank"):
		print(f"Name: {gb_record.name}, features: {len(gb_record.features)}")

#do the B. subtilis separately as the paxdb data uses the old locus tag

		if gb_record.name == "NC_000964":
		
			filename = str(gb_record.name)+".csv"
			filepath = os.path.join(my_path, "pax_db_genomes")
			complete_path = os.path.join(filepath, filename)
			file_object1 = open(complete_path, "w")
			file_object1.write("loc_tag,five_prime_seq\n")

			third_dict = {"A": 0, "T": 0, "C": 0, "G": 0}

			for feature in gb_record.features:
				if feature.type == "CDS":
					locus_tag = feature.qualifiers["old_locus_tag"][0]
					gene_sequence = feature.extract(gb_record.seq)
				
					five_prime_string = ""
					sequence = str(gene_sequence)
					n = 3
					codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]		
					five_prime = codons[1:12]
					five_prime_codon_string = "".join(five_prime)
					five_prime_string = five_prime_string + five_prime_codon_string
					
					file_object1.write(f"{locus_tag},{five_prime_string}\n")
				
					seq_3 = gene_sequence[2::3]
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
	
			file_object2.write(f"{gb_record.name},{A3},{T3},{C3},{G3},{BP},{GC3}\n")
	
			file_object1.close()
	
		else: 
		
			filename = str(gb_record.name)+".csv"
			filepath = os.path.join(my_path, "pax_db_genomes")
			complete_path = os.path.join(filepath, filename)
			file_object1 = open(complete_path, "w")
			file_object1.write("loc_tag,five_prime_seq\n")

			third_dict = {"A": 0, "T": 0, "C": 0, "G": 0}
	
			for feature in gb_record.features:
				if feature.type == "CDS":
					locus_tag = feature.qualifiers["locus_tag"][0]
					gene_sequence = feature.extract(gb_record.seq)
				
					five_prime_string = ""
					sequence = str(gene_sequence)
					n = 3
					codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]		
					five_prime = codons[1:12]
					five_prime_codon_string = "".join(five_prime)
					five_prime_string = five_prime_string + five_prime_codon_string
					
					file_object1.write(f"{locus_tag},{five_prime_string}\n")
				
					seq_3 = gene_sequence[2::3]
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
	
			file_object2.write(f"{gb_record.name},{A3},{T3},{C3},{G3},{BP},{GC3}\n")
	
			file_object1.close()
				
file_object2.close()

