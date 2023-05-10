library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in V5-prot for E. coli + VedIO 

data = read.csv("Ecoli_5prot_logodds_values.csv", h = TRUE)

expression_data = read.csv("logodds_Goodman.csv", h = TRUE)
expression_lo = expression_data[,c(1,3)]

#extract V5-prot log odds values and bind with VedIO

log_odds_data = data
NC_000913_lo = log_odds_data[,c("log_odds")]

log_odds_df = cbind(expression_lo, NC_000913_lo)

#Pearson's correlation between VedIO and V5-prot

ct = cor.test(expression_lo$log_odds, data$log_odds, method = "pearson")
rvals = ct$estimate

#PCA to calculate slope and intercept of regression line

pca = prcomp(~log_odds+NC_000913_lo, log_odds_df)
slp = with(pca, rotation[2,1] / rotation[1,1])
inte = with(pca, center[2] - slp*center[1])

#plot E. coli V5-prot against VedIO

pdf("figure5.pdf")
fig5 = ggplot(log_odds_df, aes(x=log_odds, y=NC_000913_lo, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V edIO", y = "V 5-prot") +
	geom_abline(slope = slp, intercept = inte, colour = "blue") +
	stat_cor(label.y.npc= 0.99, method = "pearson")
print(fig5)
dev.off()