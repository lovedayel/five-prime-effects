library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in Vnoise and V5-prot

data = read.csv("logodds_noise.csv", h = TRUE)

proteomics_data = read.csv("Ecoli_5prot_logodds_values.csv", h = TRUE)
proteomics_lo = proteomics_data[,c(1,3)]

log_odds_data = data
noise_lo = log_odds_data[,c("log_odds")]

#combine into one df

log_odds_df = cbind(proteomics_lo, noise_lo)

#Pearson's correlation between Vnoise and V5-prot

ct = cor.test(proteomics_lo$log_odds, data$log_odds, method = "pearson")
rvals = ct$estimate

#PCA to calculate slope and intercept of regression line

pca = prcomp(~log_odds+noise_lo, log_odds_df)
slp = with(pca, rotation[2,1] / rotation[1,1])
inte = with(pca, center[2] - slp*center[1])

#plot Vnoise against V5-prot

pdf("figure8.pdf")
fig8 = ggplot(log_odds_df, aes(x=log_odds, y=noise_lo, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V 5-prot", y = "V noise") +
	geom_abline(slope = slp, intercept = inte, colour = "blue") +
	stat_cor(method = "pearson")
print(fig8)
dev.off()
