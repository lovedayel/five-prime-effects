import pandas as pd
import math
import string

#read in csv with E. coli 5' sequences by gene
#read in excel file containing essential and non essential E. coli genes
#remove title row

seq_df = pd.read_csv("ecoli_five_prime_bygene.csv")
essential = pd.read_excel("mbo001183726st1.xlsx")
essential.columns = essential.iloc[0]
essential = essential.iloc[1:,].reindex()

#make df containing only the essential genes and merge with their 5' sequences

essential_df = essential.loc[essential["Essential"] != False ]

seq_df_dict = seq_df.set_index("Gene").transpose().to_dict(orient = "dict")
essential_df_dict = essential_df.set_index("Gene").transpose().to_dict(orient = "dict")

df_seq = pd.DataFrame(seq_df_dict)
df_essential = pd.DataFrame(essential_df_dict)

df_comb_e = pd.concat([df_seq, df_essential], axis = 0)
df_comb_e = df_comb_e.transpose()
df_comb_e = df_comb_e.dropna()
df_comb_e = df_comb_e.iloc[: , [0, 3]].copy()

#make df containing only non-essential genes and merge with their 5' sequences 

nonessential_df = essential.loc[essential["Essential"] != True ]

seq_df_dict = seq_df.set_index("Gene").transpose().to_dict(orient = "dict")
nonessential_df_dict = nonessential_df.set_index("Gene").transpose().to_dict(orient = "dict")

df_seq = pd.DataFrame(seq_df_dict)
df_nonessential = pd.DataFrame(nonessential_df_dict)

df_comb_n = pd.concat([df_seq, df_nonessential], axis = 0)
df_comb_n = df_comb_n.transpose()
df_comb_n = df_comb_n.dropna()
df_comb_n = df_comb_n.iloc[: , [0, 3]].copy() 


#log odds calculation for codons enriched in 5' ends of the essential genes

codons_essential_list = []
codons_nonessential_list = []

sequences_essential = df_comb_e["five_prime_seq"].to_list()
for seq in sequences_essential:
	n = 3
	codons_essential = [seq[i:i+n] for i in range(0, len(seq), n)]
	codons_essential_list = codons_essential_list + codons_essential
	
sequences_nonessential = df_comb_n["five_prime_seq"].to_list()
for seq in sequences_nonessential:
	n = 3
	codons_nonessential = [seq[i:i+n] for i in range(0, len(seq), n)]
	codons_nonessential_list = codons_nonessential_list + codons_nonessential
	
#open new file to output log odds ratio for each codon

file_object1 = open("logodds_essential_ecoli.csv", "w")
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

#calculate log odds ratio for the enrichment of each codon in the essential group compared to synonyms 

for aa in amino_acids:
	codons = codons_dict[aa]
	alpha = codons
	for c in alpha:
		c_essential = codons_essential_list.count(c)
		c_nonessential = codons_nonessential_list.count(c)
		not_c = [codon for codon in codons if not codon.startswith(c)]
		not_c_essential = 0
		not_c_nonessential = 0
		for nc in not_c:
			not_c_essential = not_c_essential + codons_essential_list.count(nc)
			not_c_nonessential = not_c_nonessential + codons_nonessential_list.count(nc)
	
		odds_ratio = (c_essential / c_nonessential) / (not_c_essential / not_c_nonessential)
		log_odds = round((math.log(odds_ratio)), 3)
		
		ends_with = c[2]

		file_object1.write(f"{c},{aa},{log_odds},{ends_with}\n")

file_object1.close()

