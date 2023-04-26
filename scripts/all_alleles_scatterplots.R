library(ggplot2)
library(reshape2)

lm_eqn <- function(data){
	eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
	list(
	a = format(unname(coef(m)[1]), digits = 2), 
	b = format(unname(coef(m)[2]), digits = 2),
	r2 = format(summary(m)$adj.r.squared, digits = 3)))
	as.character(as.expression(eq))
}

lm_pvalue <- function(data){
	p <- summary(m)$coefficients[,"Pr(>|t|)"][2]
	p1 <- substitute(p_value == value, list(value = format(unname(p), digits = 3)))
	as.character(as.expression(p1))
}

data = read.table("bf.mhcI.table.txt", header=T, sep="\t")

m <- lm(data$Primer2 ~ data$Primer1, data)
st = paste("r", round(cor(data$Primer2,data$Primer1, method="pearson"), digits = 3), sep=' = ')
png(file="primer1_vs_primer2.png", width=4, height=4, unit="in", res=300)
ggplot(data, aes(x=Primer1, y=Primer2)) + xlim(0,100) + ylim(0,100) + geom_point(pch=21, size=2, alpha=0.7, fill="#009E73") + theme_bw() + geom_smooth(method = "lm", colour="blue", linewidth=0.5, formula=y~x) + geom_text(x = 25, y = 98, label = lm_eqn(data), parse = TRUE, size=3, family = "Times New Roman") + geom_text(x = 20, y = 93, label = st, size=3,family = "Times New Roman") + geom_text(x = 25, y = 88, label = lm_pvalue(data), parse=TRUE, size=3,family = "Times New Roman") + theme(axis.title.x = element_text(size=12, margin=margin(20,0,0,0)), axis.title.y = element_text(size=12), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), aspect.ratio=1) + xlab("Primer1") + ylab("Primer2")
dev.off()

m <- lm(data$Primer3 ~ data$Primer1, data)
st = paste("r", round(cor(data$Primer3,data$Primer1, method="pearson"), digits = 3), sep=' = ')
png(file="primer1_vs_primer3.png", width=4, height=4, unit="in", res=300)
ggplot(data, aes(x=Primer1, y=Primer3)) + xlim(0,100) + ylim(0,100) + geom_point(pch=21, size=2, alpha=0.7, fill="#D55E00") + theme_bw() + geom_smooth(method = "lm", colour="blue", linewidth=0.5, formula=y~x) + geom_text(x = 25, y = 98, label = lm_eqn(data), parse = TRUE, size=3, family = "Times New Roman") + geom_text(x = 20, y = 93, label = st, size=3,family = "Times New Roman") + geom_text(x = 25, y = 88, label = lm_pvalue(data), parse=TRUE, size=3,family = "Times New Roman") + theme(axis.title.x = element_text(size=12, margin=margin(20,0,0,0)), axis.title.y = element_text(size=12), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), aspect.ratio=1) + xlab("Primer1") + ylab("Primer3")
dev.off()

m <- lm(data$Primer3 ~ data$Primer2, data)
st = paste("r", round(cor(data$Primer3,data$Primer2, method="pearson"), digits = 3), sep=' = ')
png(file="primer2_vs_primer3.png", width=4, height=4, unit="in", res=300)
ggplot(data, aes(x=Primer2, y=Primer3)) + xlim(0,100) + ylim(0,100) + geom_point(pch=21, size=2, alpha=0.7, fill="#0072B2") + theme_bw() + geom_smooth(method = "lm", colour="blue", linewidth=0.5, formula=y~x) + geom_text(x = 25, y = 98, label = lm_eqn(data), parse = TRUE, size=3, family = "Times New Roman") + geom_text(x = 20, y = 93, label = st, size=3,family = "Times New Roman") + geom_text(x = 25, y = 88, label = lm_pvalue(data), parse=TRUE, size=3,family = "Times New Roman") + theme(axis.title.x = element_text(size=12, margin=margin(20,0,0,0)), axis.title.y = element_text(size=12), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), aspect.ratio=1) + xlab("Primer2") + ylab("Primer3")
dev.off()
