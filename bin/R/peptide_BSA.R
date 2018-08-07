#graphs peptide BSA

library("ggplot2")

#print(getwd())

args <- commandArgs(TRUE)

BSA <- read.table(args[1])


ggplot(data = BSA, aes(x = No, y = Availability, fill = Availability))+
  geom_bar(stat = "identity")+
    theme_classic()+
  scale_x_continuous(breaks = round(seq(min(BSA$No), max(BSA$No), by = 1),1))+
  theme(legend.position="none")+
  ggtitle("Availability of peptide in bound pMHC")+
  xlab("Residue No.")
  theme(plot.title = element_text(hjust = 0.5))

ggsave("BSA.png", plot = last_plot())
write.table(BSA[,5], "peptide_BSA_piped.txt", sep="\t",
            row.names = FALSE, quote = FALSE,
            col.names = FALSE) 
