import os
from pathlib import Path
import pandas as pd
import math

#read in P. aeruginosa 5' sequence csv and P. aeruginosa protein abundance csv

my_path = os.getcwd()

path_seq = os.path.join(my_path, "pax_db_genomes/NC_002516.csv")
path_prot = os.path.join(my_path, "pax_db_proteomes/Paeruginosa_protein_abundance.csv")

seq_df = pd.read_csv(path_seq)
abundance_df = pd.read_csv(path_prot)
abundance_df = abundance_df.drop(abundance_df.columns[0], axis=1)

#combine the two CSVs into one dataframe by locus tags

seq_df_dict = seq_df.set_index("loc_tag").transpose().to_dict(orient = "dict")
abundance_df_dict = abundance_df.set_index("loc_tag").transpose().to_dict(orient = "dict")

df_seq = pd.DataFrame(seq_df_dict)
df_abundance = pd.DataFrame(abundance_df_dict)

df_comb = pd.concat([df_seq, df_abundance], axis = 0)
df_comb = df_comb.transpose()
df_comb = df_comb.dropna()

#extract top and bottom quartile of genes based on protein abundance

q75 = df_comb["abundance"].quantile(q = 0.75)
high = df_comb[df_comb["abundance"].ge(q75)]

q25 = df_comb["abundance"].quantile(q = 0.25)
low = df_comb[df_comb["abundance"].le(q25)]

#convert 5' sequences of high and low abundance genes to codons and add to their respective lists

codons_low_list = []
codons_high_list = []

sequences_low = low["five_prime_seq"].to_list()
for seq in sequences_low:
	n = 3
	codons_low = [seq[i:i+n] for i in range(0, len(seq), n)]
	codons_low_list = codons_low_list + codons_low
	
sequences_high = high["five_prime_seq"].to_list()
for seq in sequences_high:
	n = 3
	codons_high = [seq[i:i+n] for i in range(0, len(seq), n)]
	codons_high_list = codons_high_list + codons_high
	
#open new file to output log odds ratio for each codon

file_object1 = open("Paeruginosa_5prot_logodds_values.csv", "w")
file_object1.write("codon,amino_acid,log_odds,ends_with\n")

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

#calculate log odds ratio for the enrichment of each codon in the highly abundant group compared to synonyms 

for aa in amino_acids:
	codons = codons_dict[aa]
	alpha = codons
	for c in alpha:
		c_high = codons_high_list.count(c)
		c_low = codons_low_list.count(c)
		not_c = [codon for codon in codons if not codon.startswith(c)]
		not_c_high = 0
		not_c_low = 0
		for nc in not_c:
			not_c_high = not_c_high + codons_high_list.count(nc)
			not_c_low = not_c_low + codons_low_list.count(nc)
	
		odds_ratio = (c_high / c_low) / (not_c_high / not_c_low)
		log_odds = round((math.log(odds_ratio)), 3)
		
		ends_with = c[2]

		file_object1.write(f"{c},{aa},{log_odds},{ends_with}\n")

file_object1.close()