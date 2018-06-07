#graphs peptide BSA

library("ggplot2")
args <- commandArgs(TRUE)
setwd(args[1])
BSA <- read.table("peptide_BSA.txt")
BSA[,1] <- NULL
BSA[,5] <- NULL
BSA[,5] <- NULL
colnames(BSA) <- c("Amino Acid", "Residue", "ASA", "BSA")
BSA[,1] <-  substring(BSA[,1], 3)

#remember to go and check if we want 100 - score or just score
BSA$Availability <- (((BSA$ASA - BSA$BSA) / BSA$ASA) * 100)

ggplot(data = BSA, aes(x = Residue, y = Availability, fill = Availability))+
  geom_bar(stat = "identity")+
  scale_x_continuous(breaks = round(seq(min(BSA$Residue), max(BSA$Residue), by = 1),1))+
  theme(legend.position="none")+
  ggtitle("Availability of peptide in bound pMHC")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave("BSA.png", plot = last_plot())
write.table(BSA[,5], "peptide_BSA_piped.txt", sep="\t",
            row.names = FALSE, quote = FALSE,
            col.names = FALSE) 