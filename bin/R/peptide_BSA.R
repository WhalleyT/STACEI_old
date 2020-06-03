#graphs peptide BSA

library("ggplot2")
library("dplyr")

#print(getwd())

args <- commandArgs(TRUE)

BSA <- read.table(args[1])

BSA  %<>% 
  filter(Acceptor == "MHCa" | Acceptor == "MHCb") %>%
  group_by(No) %>%
  dplyr::summarize(mean_avail = mean(Availability, na.rm=TRUE))

BSA$peptide <- "peptide"

plt <- ggplot(data = BSA, aes(x = No, y = mean_avail, fill = peptide))+
  geom_bar(stat = "identity")+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(min(BSA$No), max(BSA$No), by = 1),1))+
  theme(legend.position="none")+
  ggtitle("Availability of peptide in bound pMHC")+
  xlab("Residue No.")+
  ylab("Average availability %")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("peptide" = "#FFFF00"))

ggsave("BSA.png", plot = plt)
write.table(BSA[,5], "peptide_BSA_piped.txt", sep="\t",
            row.names = FALSE, quote = FALSE,
            col.names = FALSE) 
