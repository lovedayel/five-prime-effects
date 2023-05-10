#script purpose: merge the expression log odds data with the log odds
#data for each genome into one dataframe 

expression_data = read.csv("logodds_Goodman.csv", h = TRUE)

log_odds_df = expression_data

path1 = file.path(getwd(), "genome_folders", "*", "*.csv")

files_list = Sys.glob(path1)

for (fl in files_list) {
	genome_data = read.csv(fl, h = TRUE)
	log_odds = genome_data[3]
	log_odds_df = cbind(log_odds_df, log_odds)
}

write.csv(x = log_odds_df, file = "all_log_odds.csv")

