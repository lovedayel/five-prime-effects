library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gghighlight)

#read in VedIO and VTO -> combine into one df

data = read.csv("logodds_Goodman.csv", h = TRUE)

core_data = read.csv("Ecoli_TO_logodds_values.csv", h = TRUE)

df = cbind(data, core_data)
df = df[,c(1,3,6)]

#Pearson's correlation between VedIO and VTO

ct = cor.test(df$log_odds, df$core_lo, method = "pearson")
rvals = ct$estimate

#PCA to calculate slope and intercept of regression line

pca = prcomp(~core_lo+log_odds, df)
slp = with(pca, rotation[2,1] / rotation[1,1])
inte = with(pca, center[2] - slp*center[1])

#plot VedIO against VTO

pdf("figureS1.pdf")
figS1 = ggplot(df, aes(x=core_lo, y=log_odds, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V TO", y = "V edIO") +
	geom_abline(slope = slp, intercept = inte, colour = "blue") +
	stat_cor(label.y.npc= 0.95, label.x.npc = 0.7, method = "pearson")
print(figS1)
dev.off()



