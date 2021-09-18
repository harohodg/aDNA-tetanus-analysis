library(ggplot2, lib.loc='../R/libraries')
library(ggpubr, lib.loc='../R/libraries')

tb <- read.table("mapdamage_with_metadata.csv",sep=",",header=TRUE,row.names=1)

options(repr.plot.width = 6, repr.plot.height = 4)

ggplot(tb,aes(x=MidpointCalibratedDate,y=X3pG1)) + geom_point() + geom_smooth(method=lm,color="red")+theme_bw()+ylab("3pG0")+xlab("MidpointCalibratedDate") + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.y=0.5) +  stat_regline_equation(label.y=0.6)
ggsave('sampleAge_vs_tetanusDamage_withStats_3PG0.pdf')

ggplot(tb,aes(x=MidpointCalibratedDate,y=log(X3pG1))) + geom_point() + geom_smooth(method=lm,color="red")+theme_bw()+ylab("log(3pG0)")+xlab("MidpointCalibratedDate") + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.y=0.0) +  stat_regline_equation(label.y=-0.5)
ggsave('log_sampleAge_vs_tetanusDamage_withStats_3PG0.pdf')


ggplot(tb,aes(x=MidpointCalibratedDate,y=X5pC1)) + geom_point() + geom_smooth(method=lm,color="red")+theme_bw()+ylab("5pC0")+xlab("MidpointCalibratedDate") + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.y=0.5) +  stat_regline_equation(label.y=0.6)
ggsave('sampleAge_vs_tetanusDamage_withStats_5PC0.pdf')

ggplot(tb,aes(x=MidpointCalibratedDate,y=log(X5pC1))) + geom_point() + geom_smooth(method=lm,color="red")+theme_bw()+ylab("log(5pC0)")+xlab("MidpointCalibratedDate") + stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.y=0.0) +  stat_regline_equation(label.y=-0.5)
ggsave('log_sampleAge_vs_tetanusDamage_withStats_5PC0.pdf')
