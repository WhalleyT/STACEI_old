#graphs peptide BSA
args <- commandArgs(TRUE)
library("ggplot2")
setwd(args[1])
BSA <- read.table("TCR_BSA.txt")
colnames(BSA) <- c("No.", "Chain", "Amino Acid", "Residue",
                   "ASA", "BSA", "Delta G", "Delta Gi")

alpha <- args[2]
beta <- args[3]

BSA$Availability <- (((BSA$ASA - BSA$BSA) / BSA$ASA) * 100)

#get CDR loops
fasta <- readLines(args[4])
for(i in 1:length(fasta))
{
  line <- fasta[i]
  str  <- strsplit(line, "\\|")
  if(substring(line, 1, 1) == ">")
  {
    if(str[[1]][2] == alpha)
      {
      alpha_line <- line
    }
  }
  if(substring(line, 1, 1) == ">")
  {
    if(str[[1]][2] == beta)
    {
      beta_line <- line
    }
  }
}


clean_line <- function(line)
{
  line <- strsplit(line, "\\[")[[1]][2]
  line <- gsub("\\]", "", line)
  line <- as.vector(strsplit(line, ",")[[1]])
  line <- as.numeric(line)
  
  maxline <- max(line)
  minline <- min(line)
  return(cbind(minline, maxline))
}

get_loops <- function(line)
{
  line  <- strsplit(line, "\\|")
  cdr1  <- line[[1]][7]
  cdr2  <- line[[1]][8]
  cdr3  <- line[[1]][9]
  cdrfw <- line[[1]][10]
  cdr1  <- clean_line(cdr1)
  cdr2  <- clean_line(cdr2)
  cdr3  <- clean_line(cdr3)
  cdrfw <- clean_line(cdrfw)
  return(rbind(cdr1, cdr2, cdr3, cdrfw))
}

ranges_alpha <- get_loops(alpha_line)
ranges_beta  <- get_loops(beta_line)

annotate_loops <- function(minmaxdf, df, chain, chain_str)
{
  subdf <- df[as.character(df$Chain) == chain,]
  subdf$Chain_str <- chain_str
  subdf$CDR <- NA
  
  for(i in 1:nrow(subdf))
    {
    residue <- as.numeric(as.character(subdf[i,4]))
    if(is.na(residue) == FALSE)
    {
      if(residue <= minmaxdf[1,2] && residue >= minmaxdf[1,1])
      {
        subdf[i,11] <- "CDR1"
      }
      if(residue <= minmaxdf[2,2] && residue >= minmaxdf[2,1])
      {
        subdf[i,11] <- "CDR2"
      }
      if(residue <= minmaxdf[3,2] && residue >= minmaxdf[3,1])
      {
        subdf[i,11] <- "CDR3"
      }
      if(residue <= minmaxdf[4,2] && residue >= minmaxdf[4,1])
      {
        subdf[i,11] <- "CDRfw"
      }
    }
    }
  subdf <- subdf[is.na(subdf$CDR) == FALSE,]
  subdf[is.na(subdf)] <- 0
  return(subdf)
}



alpha_plot <- annotate_loops(ranges_alpha, BSA, alpha, "alpha")
beta_plot  <- annotate_loops(ranges_beta, BSA, beta, "beta")

minus_to_start <- function(ranges, df)
{
  df$Residue <- as.numeric(as.character(df$Residue))
  df[is.na(df)] <- 0
  
  cdr1 <- df[df$CDR == "CDR1",]
  cdr1$Residue <- cdr1$Residue - (ranges[1,1] - 1)
  cdr2 <- df[df$CDR == "CDR2",]
  cdr2$Residue <- cdr2$Residue - (ranges[2,1] - 1)
  cdr3 <- df[df$CDR == "CDR3",]
  cdr3$Residue <- cdr3$Residue - (ranges[3,1] - 1)
  cdrfw <- df[df$CDR == "CDRfw",]
  cdrfw$Residue <- cdrfw$Residue - (ranges[4,1] - 1)
  out <- rbind(cdr1, cdr2, cdr3, cdrfw)
  return(out)
}

alpha_plot <- minus_to_start(ranges_alpha, alpha_plot)
beta_plot  <- minus_to_start(ranges_beta, beta_plot)

plot_df <- rbind(alpha_plot, beta_plot)
plot_df$Availability[is.na(plot_df$Availability)] <- 0
plot_df$Residue <- as.numeric(as.character(plot_df$Residue))



ggplot(data = plot_df, aes(x = Residue, y = Availability, fill = Availability))+
  geom_bar(stat = "identity")+
  theme(legend.position="none")+
  ggtitle("Availability of CDR Loops")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_grid(Chain_str ~ CDR)
ggsave("TCR_BSA_rotated.eps", plot = last_plot())

ggplot(data = plot_df, aes(x = Residue, y = Availability, fill = Availability))+
  geom_bar(stat = "identity")+
  theme(legend.position="none")+
  ggtitle("Availability of CDR Loops")+
  theme(plot.title = element_text(hjust = 0.5))+
  facet_grid(CDR ~ Chain_str)

ggsave("TCR_BSA.eps", plot = last_plot())
write.table(BSA[,5], "TCR_BSA_numbers.txt", sep="\t",
            row.names = FALSE, quote = FALSE,
            col.names = FALSE) 