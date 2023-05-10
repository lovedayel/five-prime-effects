import pandas as pd
import math
import string

#read in csv with E. coli 5' sequences by loctag
#read in excel file containing noise data for E. coli genes

seq_df = pd.read_csv("NC_000913.csv")

columns = ["B Number", "Noise_Protein"]
noise_df = pd.read_excel("TableS6.xls", usecols = columns)

#merge dataframes 

seq_df_dict = seq_df.set_index("loctag").transpose().to_dict(orient = "dict")
noise_df_dict = noise_df.set_index("B Number").transpose().to_dict(orient = "dict")

df_seq = pd.DataFrame(seq_df_dict)
df_noise = pd.DataFrame(noise_df_dict)

df_comb = pd.concat([df_seq, df_noise], axis = 0)
df_comb = df_comb.transpose()
df_comb = df_comb.dropna()

#split genes into top and bottom 40% by noise

q75 = df_comb["Noise_Protein"].quantile(q = 0.6)
high = df_comb[df_comb["Noise_Protein"].ge(q75)]

q25 = df_comb["Noise_Protein"].quantile(q = 0.4)
low = df_comb[df_comb["Noise_Protein"].le(q25)]


#log odds calculation for codons enriched in 5' ends of genes with low noise values

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

file_object1 = open("logodds_noise.csv", "w")
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

#calculate log odds ratio for the enrichment of each codon in the highly expressed group compared to synonyms 

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
	
		odds_ratio = (c_low / c_high) / (not_c_low / not_c_high)
		log_odds = round((math.log(odds_ratio)), 3)
		
		ends_with = c[2]

		file_object1.write(f"{c},{aa},{log_odds},{ends_with}\n")

file_object1.close()
