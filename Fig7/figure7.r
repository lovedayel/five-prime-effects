library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gghighlight)

#read in E. coli V5-prot and VTO

proteomics_data = read.csv("Ecoli_5prot_logodds_values.csv", h = TRUE)

core_data = read.csv("Ecoli_TO_logodds_values.csv", h = TRUE)

#combine into one df

df = cbind(proteomics_data, core_data)
df = df[,c(1,3,6)]

#Pearson's correlation between V5-prot and VTO

ct = cor.test(df$log_odds, df$core_lo, method = "pearson")
rvals = ct$estimate

#PCA to calculate slope and intercept of regression line

pca = prcomp(~core_lo+log_odds, df)
slp = with(pca, rotation[2,1] / rotation[1,1])
inte = with(pca, center[2] - slp*center[1])

#plot V5-prot against VTO

pdf("figure7.pdf")
fig7 = ggplot(df, aes(x=core_lo, y=log_odds, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V TO", y = "V 5-prot") +
	geom_abline(slope = slp, intercept = inte, colour = "blue") + 
	stat_cor(method = "pearson")
print(fig7)
dev.off()


