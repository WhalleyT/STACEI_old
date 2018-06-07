setwd("~/PycharmProjects/STACEI/coeliac")

library("ggplot2")
library("dplyr")

files <- Sys.glob("*TCR_to_pMHC_contacts_clean.txt")

# okay, so alex specifically needs ratio of CDR2a:CDR2b, but it's probably easiest to just
# do it by all CDR to pmhc and group the subsets

datalist = list()

for(f in 1:length(files))
{
  df            <- read.table(files[f], row.names = NULL, header = TRUE, sep = "\t")
  name          <- gsub("_TCR_to_pMHC_contacts_clean.txt", "", files[f])
  df$name       <- name
  datalist[[f]] <- df
}

#tidy columns
meta_contacts           <- do.call(rbind, datalist)
cols                    <- colnames(meta_contacts)
meta_contacts$row.names <- NULL
meta_contacts$Type      <- NULL
colnames(meta_contacts) <- cols[3:length(cols)]

#count
counts <- meta_contacts %>%
  group_by(name, Donor_Annotation, Acceptor_Chain, Type) %>%
  tally

#make sure it's peptide only and no non-cdrs
#counts <- counts[counts$Acceptor_Chain == "peptide",]
counts <- counts[counts$Donor_Annotation != "",]
counts <- counts[!is.na(counts$Donor_Annotation),]

#order cdrloops for plot
counts$Donor_Annotation <- factor(counts$Donor_Annotation, 
                                     levels=c("CDR1a", "CDR1b", "CDR2a", "CDR2b", 
                                              "CDR3a", "CDR3b", "FWa", "FWb"))

ggplot(counts, aes(x = `Donor_Annotation`, y = n, fill = Type))+
  geom_bar(stat = "identity")+
  ggtitle("Contact contribution of CDR loops")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title.align=  0.5)+
  scale_fill_brewer(palette = "Set2")+
  facet_grid(name ~ Acceptor_Chain)

ggplot(counts, aes(x = `Donor_Annotation`, y = n, fill = Acceptor_Chain))+
  geom_bar(stat = "identity")+
  ggtitle("Contact contribution of CDR loops")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title.align=  0.5)+
  scale_fill_brewer(palette = "Set2")+
  facet_wrap(~ name, ncol=2)

#count based on residue too
rescounts <- meta_contacts %>%
  group_by(name, Donor_Annotation, Acceptor_Chain, 
           Type, Donor_ResNum) %>%
  tally
