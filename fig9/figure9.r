library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in Vess and VedIO

data = read.csv("logodds_essential_ecoli.csv", h = TRUE)

expression_data = read.csv("logodds_Goodman.csv", h = TRUE)
expression_lo = expression_data[,c(1,3)]

log_odds_data = data
essential_lo = log_odds_data[,c("log_odds")]

#combine into one df

log_odds_df = cbind(expression_lo, essential_lo)

#Pearson's correlation between VedIO and Vess

ct = cor.test(expression_lo$log_odds, data$log_odds, method = "pearson")
rvals = ct$estimate

#PCA to calculate slope and intercept of regression line

pca = prcomp(~log_odds+essential_lo, log_odds_df)
slp = with(pca, rotation[2,1] / rotation[1,1])
inte = with(pca, center[2] - slp*center[1])

#plot Vess against VedIO

pdf("figure9.pdf")
fig9 = ggplot(log_odds_df, aes(x=log_odds, y=essential_lo, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V edIO", y = "V essential") +
	geom_abline(slope = slp, intercept = inte, colour = "blue") +
	stat_cor(method = "pearson")
print(fig9)
dev.off()