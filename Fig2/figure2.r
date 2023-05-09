library(ggplot2)
library(ggrepel)
library(ggpubr)

#read in VedIO and extract log odds values

edIO_data = read.csv("logodds_Goodman.csv", h = TRUE)
edIO_lo = edIO_data[,c(1,3)]

#read in V5-core values 

five_core_data = read.csv("log_odds_five_prime_core.csv", h = TRUE)

#extract E. coli V5-core values

NC_000913_lo = five_core_data[,c("NC_000913")]

NC_000913_df = cbind(edIO_lo, NC_000913_lo)

#PCA to calculate slope and intercept of regression line

apca = prcomp(~log_odds+NC_000913_lo, NC_000913_df)
aslp = with(apca, rotation[2,1] / rotation[1,1])
ainte = with(apca, center[2] - aslp*center[1])

#plot V5-core vs. VedIO for E. coli

pdf("figure2.pdf")
fig8 = ggplot(NC_000913_df, aes(x=log_odds, y=NC_000913_lo, label = codon)) +
	geom_point(color = "red") +
	geom_text_repel() +
	theme_bw() +
	labs(x = "V edIO", y = "V 5-core") +
	geom_abline(slope = aslp, intercept = ainte, colour = "blue") +
	stat_cor(method = "pearson")
print(fig8)
dev.off()
