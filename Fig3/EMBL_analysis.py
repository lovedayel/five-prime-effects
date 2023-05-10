import os
import re
import glob
from pathlib import Path

my_path = os.getcwd()

#make a list of all embl files 

path_embl = os.path.join(my_path, "*.embl")
files_list = glob.glob(path_embl)

#make csv file containing all accession numbers 

file_object1 = open("genomes.csv", "w")
file_object1.write("accession,bacteria,genome_size\n")

#for each genome in files_list... 
#capture accession number

for fle in files_list:
	accession = Path(fle).stem
	
#read in file
	
	genome = open(fle, "r")
	read_genome = (genome.read())
	genome.close()
	
#capture bacteria name
	
	pat1 = "XX\nOS(.*?)\n"
	bacteria = re.search(pat1, read_genome)
	bug = bacteria.group(1).strip()
	
#capture genome length
	
	top_line_search = re.search("(ID.*?)\n", read_genome)
	top_line = top_line_search.group(1)
	top_line_split = top_line.split(";")
	BP = top_line_split[-1]
	pure_BP = re.sub("[^0-9]", "", BP)
	
	file_object1.write(f"{accession},{bug},{pure_BP}\n")

#make folder named after the accession number (remove .tfa files from folder if it already exists)
#(requires a folder called 'genome_folders' to be made within the working directory)

	dir_path = os.path.join(my_path, "genome_folders", accession)
	gene_path = os.path.join(my_path, "genome_folders", accession, "*.tfa")
	try:
		os.mkdir(dir_path)
		print(f"Directory {accession} created")
	except FileExistsError:
		print(f"Directory {accession} already exists")
		existing_files = glob.glob(gene_path)
		for fl in existing_files:
			os.remove(fl)
	
#extract the sequence
	
	full_seq = read_genome.rsplit(";", 1)[1]
	
#clean the sequence (remove any white space, numbers)
	
	pure_seq = re.sub("[^a-z]", "", full_seq.lower())
	
#extract the annotation and turn into list
	
	just_annotation = read_genome.split("XX\nSQ")[0]
	pure_annotation = just_annotation.split("FT   gene")
	if len(pure_annotation) < 100:
		pure_annotation = just_annotation.split("FT   CDS")
	genelist = pure_annotation[1:]
	
#turn annotation into list
	
#for each gene in list…
#capture each gene's position and extract sequence from sequence

	gene_count = 0
	
	for gene in genelist:
		gene_count = gene_count + 1
		try:
			stripped_gene = gene.strip()
			pat = "(\d+)\.\.(\d+)"
			coordinates = re.search(pat, stripped_gene)
			unshifted_start = int(coordinates.group(1))
			shifted_start = unshifted_start -1 
			finish = int(coordinates.group(2))
			sequence = (pure_seq[shifted_start:finish])
		except:
			#print(f"error occurred with {accession} when capturing sequence)
			continue
		
#gene length?
		
		gene_length = finish - unshifted_start
		
#if complement reverse sequence & swap bases		

		complementary_bases = {"a":"t", "t":"a", "c":"g", "g":"c"}
		if stripped_gene.startswith("complement"):
			strand = "complement"
			complement_sequence = ''
			for base in sequence:
				if base in complementary_bases:
					complement_sequence = complement_sequence + complementary_bases[base]
				else:
					complement_sequence = complement_sequence + base
			sequence = complement_sequence[::-1]
		else:
			strand = "leading"
			sequence = sequence		
		
#capture locus tag 
		
		try:
			pat2 = 'locus_tag="(.+?)"'
			loctag_full = re.search(pat2, gene)
			loctag = loctag_full.group(1)
		except:
			loctag = gene_count
			
#capture translation table

		try:
			tt_pat = "/transl_table=(.*?)\n"
			tt = re.search(tt_pat, gene)
			transl_table = tt.group(1)
		except:
			transl_table = "n/a"
		
#make lists of translation table stops

		stops11 = ["taa", "tga", "tag"]
		stops4 = ["taa", "tag"]
		
#create list of genes with problems in sequence and list with which problems they contain
		
		rubbish_bin = []
		errors = []
		
#multiple of 3? 
		
		if len(sequence)%3 != 0:
			rubbish_bin.append(gene)
			errors.append("not multiple of 3")
		
#ends with a stop? 
		
		if transl_table == "11":
			if not sequence.endswith(tuple(stops11)):
				rubbish_bin.append(gene)
				errors.append("doesn't end with stop codon")
		elif transl_table == "4":
			if not sequence.endswith(tuple(stops4)):
				rubbish_bin.append(gene)
				errors.append("doesn't end with stop codon")
				
#no internal stop?
		
		n = 3
		codons = [sequence[i:i+n] for i in range(0, len(sequence), n)]
		if transl_table == "11":
			for stop in stops11:
				if stop in codons[:-1]:
					rubbish_bin.append(gene)
					errors.append("internal stop codon")
		if transl_table == "4":
			for stop in stops4:
				if stop in codons[:-1]:
					rubbish_bin.append(gene)
					errors.append("internal stop codon")
		
#starts _TG? 
		
		if sequence[1] != "t" or sequence[2] != "g":
			rubbish_bin.append(gene)
			errors.append("doesn't begin with start codon")
			
#convert error list to string
		
		errors_str = "; ".join(errors)
	
#write the sequence to a file named after locus tag and move to relevant folder
	
		if len(rubbish_bin) > 0:
			filename = str(loctag)+".tfa"
			filepath = "/Users/lovedaylewin/Documents/Bath Year 3/FYP/Python_Scripts/learning_python/filtered_genomes/rubbish/"
			complete_path = os.path.join(filepath, filename)
			with open(complete_path, "w") as f:
				f.write(f"{errors_str}\n{sequence}")
				f.close()
		else:
			filename = str(loctag)+".tfa"
			filepath = os.path.join("/Users/lovedaylewin/Documents/Bath Year 3/FYP/Python_Scripts/learning_python/filtered_genomes/genome_folders/", accession)
			complete_path = os.path.join(filepath, filename)
			with open(complete_path, "w") as f:
				f.write(f">{accession}; {bug}; {gene_length} BP; {strand}; translation table {transl_table}\n{sequence}")
				f.close()
				
file_object1.close()