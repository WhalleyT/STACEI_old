setwd("~/PycharmProjects/STACEI")

data <- read.csv("868_TCR_to_pMHC_contacts_clean.txt", sep="\t", row.names = NULL)
data <- data[data$Donor_ResNum %in% c("CDR3a", "CDR3b"),]
data <- data[data$Acceptor_Chain_Letter == "peptide",]

possible <- data %>%
  group_by(Donor_ResNum, Donor_Annotation, Donor_Chain) %>%
  tally()