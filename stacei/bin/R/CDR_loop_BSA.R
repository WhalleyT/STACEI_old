#graphs peptide BSA
args <- commandArgs(TRUE)
library("ggplot2")

setwd(args[1])

loops <- read.table("CDR_loop_BSA.txt")

ggplot(data = loops, aes(x=V1, y=V2))+
         geom_bar(stat="identity")+
  ggtitle("Total Buried Surface Area of CDR Loops")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Buried Surface Area")+
  xlab("CDR Loop")
ggsave("BSA_sum_CDR.png", plot = last_plot())