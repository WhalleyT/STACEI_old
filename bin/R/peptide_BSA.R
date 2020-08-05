#graphs peptide BSA

install_load <- function (package1, ...)  {   
  
  # convert arguments to vector
  
  # start loop to determine if each package is installed
  for(package in package1){
    
    # if package is installed locally, load
    if(package %in% rownames(installed.packages()))
      do.call('library', list(package))
    
    # if package is not installed locally, download, then load
    else {
      install.packages(package)
      do.call("library", list(package))
    }
  } 
}
packages <- c("ggplot2", "dplyr", "magrittr")

install_load(packages)

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
write.table(BSA$mean_avail, "peptide_BSA_piped.txt", sep="\t",
            row.names = FALSE, quote = FALSE,
            col.names = FALSE) 
