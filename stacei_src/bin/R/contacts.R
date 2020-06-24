#!/usr/bin/env Rscript
#this scripts generates contact graphs for TCR-pMHC
#use 1E6-RQFA as example here

#install and load dependencies
list.of.packages <- c("ggplot2", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("dplyr")
library("ggplot2")

#read in command line parameters
args <- commandArgs(trailingOnly = TRUE)
MHCa_letter      <- args[3]
MHCb_letter      <- args[4]
peptide_letter   <- args[5]
TCRa_letter      <- args[6]
TCRb_letter      <- args[7]
TCR_to_pMHC_file <- args[1]
MHC_to_pep_file  <- args[2]

setwd(args[8])
  
#read in files (change to argv later)
TCR_to_pMHC <- read.table(TCR_to_pMHC_file, sep="\t", 
                          header = TRUE, row.names = NULL)
MHC_to_pep  <- read.table(MHC_to_pep_file, sep="\t", 
                          header = TRUE, row.names = NULL)

TCR_to_pMHC$Type <- NULL
MHC_to_pep$Type <- NULL

colnames(TCR_to_pMHC) <- c(colnames(TCR_to_pMHC))[2:17]
colnames(MHC_to_pep) <- c(colnames(MHC_to_pep))[2:17]

#set booleans for what contacts actually happen
MHCa_in_MHC_to_pep     <- MHCa_letter %in%  MHC_to_pep[,1]
MHCb_in_MHC_to_pep     <- MHCb_letter %in%  MHC_to_pep[,1]
peptide_in_MHC_to_pep  <- peptide_letter %in%  MHC_to_pep[,8] 
peptide_in_TCR_to_pMHC <- peptide_letter %in% TCR_to_pMHC[,8]
TCRa_in_TCR_to_pMHC    <- TCRa_letter %in% TCR_to_pMHC[,1]
TCRb_in_TCR_to_pMHC    <- TCRb_letter %in% TCR_to_pMHC[,1]
MHCa_in_TCR_to_pMHC    <- MHCa_letter %in% TCR_to_pMHC[,8]
MHCb_in_TCR_to_pMHC    <- MHCb_letter %in% TCR_to_pMHC[,8]
######################################################################################
######################################################################################
#                                      MHC                                           #
#                                     GRAPHS                                         #
######################################################################################
######################################################################################

#MHC a contacts to peptide
if (MHCa_in_MHC_to_pep)
{
  MHCa_to_pep           <- MHC_to_pep[MHC_to_pep$Donor_Chain_Letter == MHCa_letter,]
  MHCa_to_pep_count     <- plyr::count(MHCa_to_pep, "Donor_ResNum")
  amino_acids           <- MHC_to_pep[match(MHCa_to_pep_count$Donor_ResNum, 
                                            MHCa_to_pep$Donor_ResNum),]
  amino_acids           <- dplyr::select(amino_acids, Donor_ResCode)
  MHCa_to_pep_count[,3] <- amino_acids 
  #graph it 
  ggplot(data = MHCa_to_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("MHCa to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHCa_to_pep_count$Donor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  ggsave("1.pdf")
    ggplot(data = MHCa_to_pep_count, aes(x = Donor_ResNum, y = freq))+
      geom_bar(stat="identity", aes(fill=Donor_ResCode))+
      ggtitle("MHCa to Peptide Contacts")+
      theme(plot.title = element_text(hjust = 0.5))+
      xlab("Residue Number")+
      ylab("Number of Contacts")+
      scale_x_continuous(breaks=seq(0, max(MHCa_to_pep_count$Donor_ResNum), by=5))+
      theme(axis.text.x = element_text(angle = 90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      scale_fill_discrete(name="Amino\nAcid")
    ggsave("2.pdf")
    
    MHCa_to_pep_AA           <- MHC_to_pep[MHC_to_pep$Donor_Chain_Letter == MHCa_letter,]
    MHCa_to_pep_count_AA     <- plyr::count(MHCa_to_pep, "Donor_ResCode")
    
    ggplot(MHCa_to_pep_count_AA , aes(x=1, y=freq, fill=Donor_ResCode)) +
      ggtitle("Frequency of Amino Acids Making Contact: MHCa to Peptide") +
      # black border around pie slices
      geom_bar(stat="identity", color='black') +
      # remove black diagonal line from legend
      guides(fill=guide_legend(override.aes=list(colour=NA))) +
      # polar coordinates
      coord_polar(theta='y') +
      # label aesthetics
      theme(axis.ticks=element_blank(),  # the axis ticks
            axis.title=element_blank(),  # the axis labels
            axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
            axis.text.x=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      labs(fill="Amino\nAcid")
    ggsave("3.pdf")
    
}

#MHC b contacts to peptide
if (MHCb_in_MHC_to_pep)
{
  MHCb_to_pep           <- MHC_to_pep[MHC_to_pep$Donor_Chain_Letter == MHCb_letter,]
  MHCb_to_pep_count     <- plyr::count(MHCb_to_pep, "Donor_ResNum")
  amino_acids           <- MHC_to_pep[match(MHCb_to_pep_count$Donor_ResNum, 
                                            MHCb_to_pep$Donor_ResNum),]
  amino_acids           <- dplyr::select(amino_acids, Donor_ResCode)
  MHCb_to_pep_count[,3] <- amino_acids 
  #graph it 
  ggplot(data = MHCb_to_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("MHCb to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHCb_to_pep_count$Donor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("4.pdf")
  
  ggplot(data = MHCb_to_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=Donor_ResCode))+
    ggtitle("MHCb to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHCb_to_pep_count$Donor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_discrete(name="Amino\nAcid")
  ggsave("5.pdf")
  
  MHCb_to_pep_AA           <- MHC_to_pep[MHC_to_pep$Donor_Chain_Letter == MHCb_letter,]
  MHCb_to_pep_count_AA     <- plyr::count(MHCb_to_pep, "Donor_ResCode")
  ggplot(MHCb_to_pep_count_AA , aes(x=1, y=freq, fill=Donor_ResCode)) +
    ggtitle("Frequency of Amino Acids Making Contact: MHCb to Peptide") +
    # black border around pie slices
    geom_bar(stat="identity", color='black') +
    # remove black diagonal line from legend
    guides(fill=guide_legend(override.aes=list(colour=NA))) +
    # polar coordinates
    coord_polar(theta='y') +
    # label aesthetics
    theme(axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
          axis.text.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    labs(fill="Amino\nAcid")
  ggsave("6.pdf")
  
}

if (MHCa_in_MHC_to_pep && MHCb_in_MHC_to_pep)
{
  MHCa_to_pep_count[,4] <- "alpha"
  MHCb_to_pep_count[,4] <- "beta"
  MHC_contacts_to_pep <- rbind(MHCa_to_pep_count, MHCb_to_pep_count)
  colnames(MHC_contacts_to_pep)[4] <- "MHC chain"
  
  ggplot(data = MHC_contacts_to_pep, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    facet_grid(`MHC chain` ~ .)+
    ggtitle("MHC to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHC_contacts_to_pep$Donor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("7.pdf")
}
  
######################################################################################
######################################################################################
#                                     PEPTIDE                                        #
#                                     GRAPHS                                         #
######################################################################################
######################################################################################

#ratio of MHC:TCR contacts
if(peptide_in_MHC_to_pep && peptide_in_TCR_to_pMHC)
{
  Pep_MHC   <- MHC_to_pep[MHC_to_pep$Acceptor_Chain_Letter == peptide_letter,]
  Pep_TCR   <- TCR_to_pMHC[TCR_to_pMHC$Acceptor_Chain_Letter == peptide_letter,]
  
  no_pep_MHC_contacts <- as.numeric(nrow(Pep_MHC))
  no_pep_TCR_contacts <- as.numeric(nrow(Pep_TCR))
  
  TCRMHC_pep    <- c("TCR", "MHC")
  counts        <- c(no_pep_TCR_contacts, no_pep_MHC_contacts)
  TCRMHC_counts <- as.data.frame(cbind(TCRMHC_pep, counts))
  
  ggplot(TCRMHC_counts , aes(x=1, y=counts, fill=TCRMHC_pep)) +
    ggtitle("Ratio of Contacts to Peptide") +
    # black border around pie slices
    geom_bar(stat="identity", color='black') +
    # remove black diagonal line from legend
    guides(fill=guide_legend(override.aes=list(colour=NA))) +
    # polar coordinates
    coord_polar(theta='y') +
    # label aesthetics
    theme(axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
          axis.text.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_fill_discrete(name="Complex")+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave("8.pdf")
  
  #master plot of peptide to whole complex
  pep_to_MHC_count     <- plyr::count(Pep_MHC, "Acceptor_ResNum")
  pep_to_TCR_count     <- plyr::count(Pep_TCR, "Acceptor_ResNum")
  pep_to_MHC_count[,3] <- "MHC"
  pep_to_TCR_count[,3] <- "TCR"
  
  
  colnames(pep_to_TCR_count)[1] <- "Acceptor_ResNum"
  peptide_contacts     <- rbind(pep_to_MHC_count, pep_to_TCR_count)
  
  ggplot(data = peptide_contacts, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("Peptide to TCR & MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(peptide_contacts$Acceptor_ResNum), by=1))+
    scale_fill_discrete(name="Complex")
  ggsave("9.pdf")
  
  ggplot(data = peptide_contacts, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    facet_grid(V3 ~ .)+
    ggtitle("Peptide to TCR & MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(peptide_contacts$Acceptor_ResNum), by=1))+
    scale_fill_discrete(name="Complex")
  ggsave("10.pdf")
}

#TCRa contacts
if(peptide_in_TCR_to_pMHC && TCRa_in_TCR_to_pMHC )
{
  TCRa_cont  <- TCR_to_pMHC[TCR_to_pMHC$Donor_Chain_Letter == TCRa_letter,]
  pep_TCRa   <- TCRa_cont[TCRa_cont$Acceptor_Chain_Letter == peptide_letter,]
  pep_TCRa_count <- plyr::count(pep_TCRa, "Acceptor_ResNum")
  TCRa_pep_count <- plyr::count(pep_TCRa, "Donor_ResNum")
  
  
  pep_TCRa_count_AA <- pep_TCRa[match(pep_TCRa_count$Acceptor_ResNum,
                                      pep_TCRa$Acceptor_ResNum),]
  
  TCRa_pep_count_AA <- pep_TCRa[match(TCRa_pep_count$Donor_ResNum, pep_TCRa$Donor_ResNum),] 
  TCRa_pep_count_AA <- TCRa_pep_count_AA[,5] 
  pep_TCRa_count_AA <- pep_TCRa_count_AA[,12] 

  pep_TCRa_count[,3] <- pep_TCRa_count_AA
  TCRa_pep_count[,3] <- TCRa_pep_count_AA
  
  TCRa_pep_count[,4] <- NA

  
  ggplot(data = pep_TCRa_count, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("Peptide to TCRa Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(pep_TCRa_count$Acceptor_ResNum), by=1))
  ggsave("11.pdf")
  
  ggplot(data = pep_TCRa_count, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("Peptide to TCRa Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(pep_TCRa_count$Acceptor_ResNum), by=1))+
    scale_fill_discrete(name="Amino Acid")
  ggsave("12.pdf")
  
  ggplot(data = TCRa_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("TCRa to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRa_pep_count$Donor_ResNum), by=5))
  ggsave("13.pdf")
  
  ggplot(data = TCRa_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("TCRa to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRa_pep_count$Donor_ResNum), by=5))+
    scale_fill_discrete(name="Amino Acid")
  ggsave("14.pdf")

}


if(peptide_in_TCR_to_pMHC &&TCRb_in_TCR_to_pMHC)
{
  TCRb_cont  <- TCR_to_pMHC[TCR_to_pMHC$Donor_Chain_Letter == TCRb_letter,]
  pep_TCRb   <- TCRb_cont[TCRb_cont$Acceptor_Chain_Letter == peptide_letter,]
  pep_TCRb_count <- plyr::count(pep_TCRb, "Acceptor_ResNum")
  TCRb_pep_count <- plyr::count(pep_TCRb, "Donor_ResNum")
  
  pep_TCRb_count_AA <- pep_TCRb[match(pep_TCRb_count$Acceptor_ResNum,
                                      pep_TCRb$Acceptor_ResNum),]
  TCRb_pep_count_AA <- pep_TCRb[match(TCRb_pep_count$Donor_ResNum, pep_TCRb$Donor_ResNum),] 
  TCRb_pep_count_AA <- TCRb_pep_count_AA[,5] 
  pep_TCRb_count_AA <- pep_TCRb_count_AA[,12] 
  
  pep_TCRb_count[,3] <- pep_TCRb_count_AA
  TCRb_pep_count[,3] <- TCRb_pep_count_AA
  
  ggplot(data = pep_TCRb_count, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("Peptide to TCRb Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(pep_TCRb_count$Acceptor_ResNum), by=1))
  ggsave("15.pdf")
  
  ggplot(data = pep_TCRb_count, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("Peptide to TCRb Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(pep_TCRb_count$Acceptor_ResNum), by=1))+
    scale_fill_discrete(name="Amino Acid")
  ggsave("16.pdf")
  
  ggplot(data = TCRb_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("TCRb to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")#+
    #scale_x_continuous(breaks=seq(0, max(TCRb_pep_count$ResNum), by=5))
  ggsave("17.pdf")
  
  ggplot(data = TCRb_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("TCRb to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    #scale_x_continuous(breaks=seq(0, max(TCRb_pep_count$ResNum), by=10))+
    scale_fill_discrete(name="Amino Acid")
  ggsave("18.pdf")
}

if(exists("TCRa_pep_count"))
{
  for(i in 1:nrow(TCRa_pep_count))
  {
    if(TCRa_pep_count[i,1] >= 23 && TCRa_pep_count[i,1] <= 33 )
    {
      TCRa_pep_count[i,4] <- "CDR1a"
    }
    else if(TCRa_pep_count[i,1] >= 48 && TCRa_pep_count[i,1] <= 58 )
    {
      TCRa_pep_count[i,4] <- "CDR2a"      
    }
    else if(TCRa_pep_count[i,1] >= 90 && TCRa_pep_count[i,1] <= 110 )
    {
      TCRa_pep_count[i,4] <- "CDR3a"       
    }    
    else if(TCRa_pep_count[i,1] >= 64 && TCRa_pep_count[i,1] <= 74 )
    {
      TCRa_pep_count[i,4] <- "CDR3fwa"         
    }
    else
    {
      TCRa_pep_count[i,4] <- "none"   
    }
  }
  
  ggplot(data = TCRa_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V4))+
    ggtitle("TCRa to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRa_pep_count$Donor_ResNum), by=5))+
    scale_fill_discrete(name="CDR Loop")
  ggsave("19.pdf")
}

if(exists("TCRb_pep_count"))
{
  for(i in 1:nrow(TCRb_pep_count))
  {
    if(TCRb_pep_count[i,1] >= 23 && TCRb_pep_count[i,1] <= 33 )
      {
      TCRb_pep_count[i,4] <- "CDR1b"
      }
    else if(TCRb_pep_count[i,1] >= 46 && TCRb_pep_count[i,1] <= 56 )
      {
      TCRb_pep_count[i,4] <- "CDR2b"      
      }
    else if(TCRb_pep_count[i,1] >= 90 && TCRb_pep_count[i,1] <= 110 )
      {
      TCRb_pep_count[i,4] <- "CDR3b"       
      }    
    else if(TCRb_pep_count[i,1] >= 64 && TCRb_pep_count[i,1] <= 74 )
      {
      TCRb_pep_count[i,4] <- "CDR3fwb"         
      }
    else
      {
      TCRb_pep_count[i,4] <- "none"   
      }
  }
  
  ggplot(data = TCRb_pep_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V4))+
    ggtitle("TCRb to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    #scale_x_continuous(breaks=seq(0, max(TCRb_pep_count$ResNum), by=10))+
    scale_fill_discrete(name="CDR Loop")
  ggsave("20.pdf")
  }


if(exists("TCRb_pep_count") && exists("TCRa_pep_count"))
{
  TCRb_pep_count[,5] <- "TCRb"
  TCRa_pep_count[,5] <- "TCRa"
  colnames(TCRb_pep_count) <- colnames(TCRa_pep_count)
  TCR_pep_contacts <- rbind(TCRa_pep_count, TCRb_pep_count)
  
  ggplot(data = TCR_pep_contacts, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V4))+
    facet_grid(V5 ~ .)+
    ggtitle("TCR to Peptide Contacts")+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCR_pep_contacts$Donor_ResNum), by=5))+
    scale_fill_discrete(name="CDR Loop")
  ggsave("21.pdf")
}

#TCR to MHC and vice versa
if(TCRa_in_TCR_to_pMHC)
{
  TCRa_to_MHC <- TCRa_cont[TCRa_cont$Acceptor_Chain_Letter ==  MHCa_letter,]
  TCRa_to_MHC_count <- plyr::count(TCRa_to_MHC, "Donor_ResNum")
  MHC_to_TCRa_count <- plyr::count(TCRa_to_MHC, "Acceptor_ResNum")
  
  amino_acids_TCR   <- TCRa_to_MHC[match(TCRa_to_MHC_count$Donor_ResNum,
                                         TCRa_to_MHC$Donor_ResNum),]
  
  amino_acids_MHC   <- TCRa_to_MHC[match(MHC_to_TCRa_count$Acceptor_ResNum,
                                         TCRa_to_MHC$Acceptor_ResNum),]
  
  amino_acids_MHC       <- amino_acids_MHC[,5]
  amino_acids_TCR       <-  amino_acids_TCR[,12]
  TCRa_to_MHC_count[,3] <- amino_acids_TCR
  MHC_to_TCRa_count[,3] <- amino_acids_MHC
  
  ggplot(data = TCRa_to_MHC_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("TCRa to MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRa_to_MHC_count$Donor_ResNum), by=5))+
    scale_fill_discrete(name="Amino Acid")+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("22.pdf")
  
  ggplot(data = TCRa_to_MHC_count, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("TCRa to MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRa_to_MHC_count$Donor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("23.pdf")
  
  ggplot(data = MHC_to_TCRa_count, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("MHC to TCRa Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHC_to_TCRa_count$Acceptor_ResNum), by=5))+
    scale_fill_discrete(name="Amino Acid")+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()))
  ggsave("24.pdf")
  
  ggplot(data = MHC_to_TCRa_count, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("MHC to TCRa Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHC_to_TCRa_count$Acceptor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()))
  ggsave("25.pdf")
}

if(TCRb_in_TCR_to_pMHC)
{
  TCRb_to_MHC <- TCRb_cont[TCRb_cont$Chain ==  MHCa_letter,]
  TCRb_to_MHC_count <- plyr::count(TCRb_to_MHC, "Donor_ResNum")
  MHC_to_TCRb_count <- plyr::count(TCRb_to_MHC, "Acceptor_ResNum")
  
  amino_acids_TCR   <- TCRb_to_MHC[match(TCRb_to_MHC_count$Donor_ResNum,
                                         TCRb_to_MHC$Donor_ResNum),]
  
  amino_acids_MHC   <- TCRb_to_MHC[match(MHC_to_TCRb_count$Acceptor_ResNum,
                                         TCRb_to_MHC$Acceptor_ResNum),]
  
  amino_acids_MHC       <- amino_acids_MHC[,5]
  amino_acids_TCR       <-  amino_acids_TCR[,12]
  TCRb_to_MHC_count[,3] <- amino_acids_TCR
  MHC_to_TCRb_count[,3] <- amino_acids_MHC
  
  ggplot(data = TCRb_to_MHC_count, aes(x = ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("TCRb to MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRb_to_MHC_count$ResNum), by=5))+
    scale_fill_discrete(name="Amino Acid")+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("26.pdf")
  
  ggplot(data = TCRb_to_MHC_count, aes(x = ResNum, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("TCRb to MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCRb_to_MHC_count$ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("27.pdf")
  
  ggplot(data = MHC_to_TCRb_count, aes(x = ResNum.1, y = freq))+
    geom_bar(stat="identity", aes(fill=V3))+
    ggtitle("MHC to TCRb Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHC_to_TCRb_count$ResNum.1), by=5))+
    scale_fill_discrete(name="Amino Acid")+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("28.pdf")
  
  ggplot(data = MHC_to_TCRb_count, aes(x = ResNum.1, y = freq))+
    geom_bar(stat="identity")+
    ggtitle("MHC to TCRb Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHC_to_TCRb_count$ResNum.1), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  ggsave("29.pdf")
}

if(TCRa_in_TCR_to_pMHC && TCRb_in_TCR_to_pMHC)
{
  MHC_to_TCRb_count[,4] <- "TCRb"
  TCRb_to_MHC_count[,4] <- "TCRb"
  
  MHC_to_TCRa_count[,4] <- "TCRa"
  TCRa_to_MHC_count[,4] <- "TCRa"
  
  colnames(MHC_to_TCRb_count) <- colnames(MHC_to_TCRa_count)
  colnames(TCRb_to_MHC_count) <- colnames(TCRa_to_MHC_count)
  
  MHC_to_TCR <- rbind(MHC_to_TCRa_count, MHC_to_TCRb_count)
  TCR_to_MHC <- rbind(TCRa_to_MHC_count, TCRb_to_MHC_count)
  
  ggplot(data = MHC_to_TCR, aes(x = Acceptor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V4))+
    ggtitle("MHC to TCR Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(MHC_to_TCR$Acceptor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_fill_discrete(name="Acceptor_Chain")
  ggsave("30.pdf")
  
  ggplot(data = TCR_to_MHC, aes(x = Donor_ResNum, y = freq))+
    geom_bar(stat="identity", aes(fill=V4))+
    ggtitle("TCR to MHC Contacts")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Residue Number")+
    ylab("Number of Contacts")+
    scale_x_continuous(breaks=seq(0, max(TCR_to_MHC$Donor_ResNum), by=5))+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_fill_discrete(name="Donor_Chain")
  ggsave("31.pdf")
}

#total CDR bars to peptide
if(peptide_in_TCR_to_pMHC)
  {
    TCR_to_pep <- TCR_to_pMHC[TCR_to_pMHC$Acceptor_Chain_Letter == peptide_letter,]
    TCR_to_pep[,17] <- NA
    for(i in 1:nrow(TCR_to_pep))
    {
      if(TCR_to_pep[i,3] >= 23 && TCR_to_pep[i,3] <= 33 )
      {
        TCR_to_pep[i,17] <- "CDR1a"
      }
      else if(TCR_to_pep[i,3] >= 48 && TCR_to_pep[i,3] <= 58 )
      {
        TCR_to_pep[i,17] <- "CDR2a"      
      }
      else if(TCR_to_pep[i,3] >= 90 && TCR_to_pep[i,3] <= 110 )
      {
        TCR_to_pep[i,17] <- "CDR3a"       
      }    
      else if(TCR_to_pep[i,3] >= 64 && TCR_to_pep[i, 3] <= 74 )
      {
        TCR_to_pep[i,17] <- "CDR3fwa"         
      }
      else
      {
        TCR_to_pep[i,17] <- "none"   
      }
    }
    
    CDR_to_pep_count <- plyr::count(TCR_to_pep, "V17")
    ggplot(CDR_to_pep_count , aes(x=1, y=freq, fill=V17)) +
      ggtitle("Ratio of CDR to Peptide Contacts") +
      # black border around pie slices
      geom_bar(stat="identity", color='black') +
      # remove black diagonal line from legend
      guides(fill=guide_legend(override.aes=list(colour=NA))) +
      # polar coordinates
      coord_polar(theta='y') +
      # label aesthetics
      theme(axis.ticks=element_blank(),  # the axis ticks
            axis.title=element_blank(),  # the axis labels
            axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
            axis.text.x=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())+
      scale_fill_discrete(name="CDR Loop")+
      theme(plot.title = element_text(hjust = 0.5))
    ggsave("32.pdf")
  }

#total CDR bars to MHC
if(MHCa_in_TCR_to_pMHC)
{
  TCR_to_pep <- TCR_to_pMHC[TCR_to_pMHC$Acceptor_Chain_Letter != peptide_letter,]
  TCR_to_pep[,17] <- NA
  for(i in 1:nrow(TCR_to_pep))
  {
    if(TCR_to_pep[i,3] >= 23 && TCR_to_pep[i,3] <= 33 )
    {
      TCR_to_pep[i,17] <- "CDR1a"
    }
    else if(TCR_to_pep[i,3] >= 48 && TCR_to_pep[i,3] <= 58 )
    {
      TCR_to_pep[i,17] <- "CDR2a"      
    }
    else if(TCR_to_pep[i,3] >= 90 && TCR_to_pep[i,3] <= 110 )
    {
      TCR_to_pep[i,17] <- "CDR3a"       
    }    
    else if(TCR_to_pep[i,3] >= 64 && TCR_to_pep[i,3] <= 74 )
    {
      TCR_to_pep[i,17] <- "CDR3fwa"         
    }
    else
    {
      TCR_to_pep[i,17] <- "none"   
    }
  }
  
  CDR_to_pep_count <- plyr::count(TCR_to_pep, "V17")
  ggplot(CDR_to_pep_count , aes(x=1, y=freq, fill=V17)) +
    ggtitle("Ratio of CDR to MHC Contacts") +
    # black border around pie slices
    geom_bar(stat="identity", color='black') +
    # remove black diagonal line from legend
    guides(fill=guide_legend(override.aes=list(colour=NA))) +
    # polar coordinates
    coord_polar(theta='y') +
    # label aesthetics
    theme(axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
          axis.text.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_fill_discrete(name="CDR Loop")+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave("33.pdf")
}
#ratio of each CDR loop to pMHC

  TCR_to_pep <- TCR_to_pMHC
  TCR_to_pep[,17] <- NA
  for(i in 1:nrow(TCR_to_pep))
  {
    if(TCR_to_pep[i,3] >= 23 && TCR_to_pep[i,3] <= 33 )
    {
      TCR_to_pep[i,17] <- "CDR1a"
    }
    else if(TCR_to_pep[i,3] >= 48 && TCR_to_pep[i,3] <= 58 )
    {
      TCR_to_pep[i,17] <- "CDR2a"      
    }
    else if(TCR_to_pep[i,3] >= 90 && TCR_to_pep[i,3] <= 110 )
    {
      TCR_to_pep[i,17] <- "CDR3a"       
    }    
    else if(TCR_to_pep[i,3] >= 64 && TCR_to_pep[i,3] <= 74 )
    {
      TCR_to_pep[i,17] <- "CDR3fwa"         
    }
    else
    {
      TCR_to_pep[i,17] <- "none"   
    }
  }
  
  CDR_to_pep_count <- plyr::count(TCR_to_pep, "V17")
  ggplot(CDR_to_pep_count , aes(x=1, y=freq, fill=V17)) +
    ggtitle("Ratio of CDR to all Contacts") +
    # black border around pie slices
    geom_bar(stat="identity", color='black') +
    # remove black diagonal line from legend
    guides(fill=guide_legend(override.aes=list(colour=NA))) +
    # polar coordinates
    coord_polar(theta='y') +
    # label aesthetics
    theme(axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
          axis.text.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_fill_discrete(name="CDR Loop")+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave("34.pdf")