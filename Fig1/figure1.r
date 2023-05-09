library(ggplot2)

#read in VedIO

data = read.csv("logodds_Goodman.csv", h = TRUE)

#plot VedIO by codon, coloured based on 3rd position nucleotide

pdf("figure1.pdf", width = 10, height = 5)
fig1 = ggplot(data, aes(x = codon, y = log_odds)) + 
	geom_col(aes(fill = ends_with)) + 
	facet_grid(~amino_acid, scales = "free_x", space = "free_x") + 
	scale_x_discrete(guide = guide_axis(angle = -90)) + 
	labs(x = "Codon", y = "Odds ratio", fill = "3rd Position \nNucleotide") + 
	theme(strip.text.x = element_text(face = "bold"), strip.background = element_rect(color = "black")) 
print(fig1)
dev.off()

