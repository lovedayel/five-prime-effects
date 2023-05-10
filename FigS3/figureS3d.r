library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in VedIO and P. aeruginosa V5-prot

data = read.csv("Paeruginosa_5prot_logodds_values.csv", h = TRUE)

expression_data = read.csv("logodds_Goodman.csv", h = TRUE)
expression_lo = expression_data[,c(1,3)]

log_odds_data = data
Paeruginosa_lo = log_odds_data[,c("log_odds")]

#combine in one df

log_odds_df = cbind(expression_lo, Paeruginosa_lo)

#Pearson's correlation between VedIO and V5-prot

ct = cor.test(expression_lo$log_odds, data$log_odds, method = "pearson")
rvals = ct$estimate

#PCA to calculate slope and intercept of regression line

pca = prcomp(~log_odds+Paeruginosa_lo, log_odds_df)
slp = with(pca, rotation[2,1] / rotation[1,1])
inte = with(pca, center[2] - slp*center[1])

#plot V5-prot against VedIO

pdf("figureS3d.pdf")
figS3d = ggplot(log_odds_df, aes(x=log_odds, y=Paeruginosa_lo, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V edIO", y = "V 5-prot") +
	geom_abline(slope = slp, intercept = inte, colour = "blue") +
	stat_cor(method = "pearson")
print(figS3d)
dev.off()

