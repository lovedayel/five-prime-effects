import os
import glob
import re
from pathlib import Path
import pandas as pd
from io import StringIO

#make a list of the PaxDb abundance data files for the five genomes

my_path = os.getcwd()

path_txt = os.path.join(my_path, "pax_db_proteomes", "*.txt")
files_list = glob.glob(path_txt)

for fle in files_list:

	name = Path(fle).stem
	
#create file path for the output of each genome
	
	filename = name+".csv"
	filepath = os.path.join(my_path, "pax_db_proteomes")
	complete_path = os.path.join(filepath, filename)

	proteome = open(fle, "r")
	read_proteome = (proteome.read())
	proteome.close()

#extract locus tag and abundance values and input into a dataframe 

	data = read_proteome.rsplit("#", 1)[1]

	df = pd.read_csv(StringIO(data), sep='\t')
	
	loctag = df["string_external_id"].str.split(".", expand = True)[1]
	
	abundance = df["abundance"].astype(int)
	
	df2 = pd.concat([loctag, abundance], axis = 1)
	df2.columns = ["loc_tag", "abundance"]
	print(df2)

#output dataframe to csv for each genome
	
	df2.to_csv (complete_path, header = True) 

