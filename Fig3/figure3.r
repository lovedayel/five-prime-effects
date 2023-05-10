#read in csv with VedIO + V5-core for all 650 bacteria

log_odds_data = read.csv("all_log_odds.csv", h = TRUE)

#extract column with VedIO

exp.lo = log_odds_data$log_odds

ncols = dim(log_odds_data)[2]

r.vals = c()
p.vals = c()

#calculate r and p values for each genome's V5-core correlated against VedIO

for (i in c(5:ncols)) {
	lo = unlist(log_odds_data[i])
	ct = cor.test(exp.lo, lo, method = "pearson")
	ct.r = ct$estimate
	ct.p = ct$p.value
	r.vals = c(r.vals, ct.r)
	p.vals = c(p.vals, ct.p)
}

#find the points that demarcate statistical significance 

accession = colnames(log_odds_data)[5:ncols]

df.cor = data.frame(acc = accession, rvals = r.vals, pvals = p.vals)
#includes bonferroni correction for multiple testing
non_sig = which(df.cor$pvals > (0.05/ncols))
df.non_sig = df.cor$rvals[non_sig]
max.r = max(df.non_sig)
min.r = min(df.non_sig)

#make output csv

write.csv(x = df.cor, file = "correlation_df.csv")

#plot histogram of R values 

pdf("figure3.pdf")
hist(r.vals,
main = "",
xlab = "Pearson's r",
col = "blue",
breaks = 50)

#add vertical lines to show which genomes show significant correlation (both + and -)

abline(v = max.r, col = "red")
abline(v = min.r, col = "red")

dev.off()

#what proportion of genomes show significant correlation?

tot_non_sig = length(non_sig)
tot = dim(df.cor)[1]
prop_sig = (tot - tot_non_sig) / tot

#which genomes show significant + correlation and which significant - correlation?

sig_pos = df.cor$acc[df.cor$rvals > 0 & df.cor$pvals < (0.05/ncols)]
sig_neg = df.cor$acc[df.cor$rvals < 0 & df.cor$pvals < (0.05/ncols)]



