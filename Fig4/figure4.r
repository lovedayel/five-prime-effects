library(ggplot2)

#read in CSV containing Pearson's r + GC3 content for each bacteria

data = read.csv("rval_vs_GC3.csv", h = TRUE)

#add column containing GC3^2 for quadratic model of r vals

data2 = transform(data, GC3_squared=GC3_content^2)

quadratic_model = lm(data2$rvals ~ data2$GC3_content + data2$GC3_squared)

#extract E. coli data point so it can be highlighted in the figure

ecoli = subset(data2, Accession == "AP010953")

#plot r values against GC3 content including quadratic model

pdf("figure4.pdf", width = 8, height = 6)
fig4 = ggplot(data2, aes(x = GC3_content, y = rvals)) +
	geom_point(colour = "blue", pch = 21) + 
	geom_point(data = ecoli, colour = "red", pch = 18, cex = 4) +
	geom_text(data = ecoli, label = "E. coli", nudge_x = 5, nudge_y = 0.05, colour = "red") +
	theme_bw() + 
	geom_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "black", se = FALSE, linewidth = 0.7) +
	labs(x = "GC3 content (%)", y = "Pearson's r'") + 
	annotate("text", x=29, y=0.85, label= "Adjusted R-squared: 0.6569, p-value: < 2.2e-16")
print(fig4)
dev.off()

