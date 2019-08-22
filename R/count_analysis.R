library(gridExtra)
library(grid)
require(phyloseq)
library(stringr)
library(ade4)
library(ggplot2)
library(reshape2)
library(reshape)
library(plyr)
library(dplyr)
library(vegan)

figsize = function(width, height) options(repr.plot.width=width, repr.plot.height=height)

KCAPplot2<- function(variable, sample, T, P, legend='variable', c.dist = c(0, 0.5, 1), draw=TRUE, ylim1=c(0,2), ylim2=c(0,2)){
df = data.frame(variable=variable, sample=sample, T=T, P=P)
df1 <- ddply(df, .variables = c("sample"), summarise, variable_m = mean(variable))
df2 <- merge(df1, df)

for_plot <- ddply(df2, .variable=c("T", "P"), summarise, avg=mean(variable/variable_m))

for_plot_P <- ddply(df2, .variable=c("P"), summarise, 
                    avg=mean(variable/variable_m), se = sd(variable/variable_m)/sqrt(length(variable)))
for_plot_T <- ddply(df2, .variable=c("T"), summarise, 
                    avg=mean(variable/variable_m), se = sd(variable/variable_m)/sqrt(length(variable)))

#for_plotf

hm <- ggplot(for_plot, aes(x=T, y=P)) + geom_tile(aes(fill=avg)) + 
theme_classic() + scale_fill_gradientn(colors = c("green", "yellow", "red"), 
                                       values= c.dist)+
theme(legend.position = "bottom", legend.title = element_blank(),
      legend.key.height = unit(2, "cm"),legend.key.width=unit(1.5, "cm"))+ 
theme(plot.margin=margin(2,2,3,2, unit="cm")) + 
ylab("percentage of maximum difference") + xlab("trimming threshold")

bar3 <- ggplot(for_plot_T, aes(x=factor(T),# fill=avg,
                           y=avg)) + 
geom_errorbar(aes(ymin=avg-se,ymax=avg+se), position=position_dodge(), width=0.4)+
geom_point()+ 
    #geom_bar(stat="identity")+
theme_classic() + ylab(paste("normalized", legend)) +
#annotate("text", label="*", x=1, y=0.998, size=7)+
xlab("trimming threshold")+ 
theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(ylim =ylim1)

bar4 <- ggplot(for_plot_P, aes(x=factor(P), 
                           #fill=avg, 
                           y=avg)) + 
geom_errorbar(aes(ymin=avg-se,ymax=avg+se), position=position_dodge(), width=0.4)+
#geom_bar(stat="identity")+
geom_point()+ 
theme_classic() + ylab(paste("normalized", legend)) +
#    coord_cartesian(ylim =ylim2) +
xlab("percentage of mismatch") + coord_flip(ylim=ylim2) +#ylim=c(0.575,1.535)) +  
#scale_fill_gradientn(colors = c("yellow","red"), values = c(0, 1))+
theme(plot.title = element_text(hjust = 0.5)) 
tmp <- ggplot_gtable(ggplot_build(hm))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]

g1 <- ggplotGrob(bar3+ theme(legend.position = "none")+ xlab(""))
g2 <- ggplotGrob(bar4 + theme(legend.position = "none")+ xlab(""))
g3 <- ggplotGrob(hm+theme(legend.position = "none"))
#g4 <- ggplotGrob(legend)

g1$widths -> g2$widths
g1$widths -> g3$widths
g2$heights -> g1$heights
g3$heights <- g1$heights

g <- arrangeGrob(g1, legend, g3, g2, ncol=2)
if (draw) grid.draw(g)
return (g)
                }

figsize(10, 8)

chi_count = read.csv("chi_count_gh.csv", row.names=NULL)

chi_count$T = str_pad(chi_count$T, 2, pad="0" )
chi_count$P = str_pad(chi_count$P, 2, pad="0")

chi_count_real <- chi_count[chi_count$template=="real",]
chi_count_template <- chi_count[chi_count$template=="template",]
chi_count_perfect <- chi_count[chi_count$template=="perfect",]

figsize(6, 5)
g <- KCAPplot2(variable = chi_count_real$count, sample = chi_count_real$Sample, 
               T=chi_count_real$T, P=chi_count_real$P, legend = "count", c.dist = c(0, 0.5,1), 
               ylim1=c(0.75, 1.09), ylim2=(c(0.25, 1.15)))

print(g)

