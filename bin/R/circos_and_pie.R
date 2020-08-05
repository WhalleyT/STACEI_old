options(show.error.locations = TRUE)

packages <- c("dplyr", "tidyr", "ggplot2", "purrr", "RColorBrewer",
              "circlize", "reshape2")


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

install_load(packages)

make_pallete <- function(df)
{
  lightblue    <- rgb(0.75, 0.75, 1.0)
  palegreen    <- rgb(0.65, 0.9, 0.65)
  yellow       <- rgb(1.0 , 1.0,	0.0)
  grey60       <- rgb(0.6 , 0.6,	0.6)
  grey80       <- rgb(0.8 , 0.8,	0.8)
  tv_red       <- rgb(1.0, 0.2, 0.2)
  tv_green     <- rgb(0.2, 1.0, 0.2)
  purpleblue   <- rgb(0.5, 0.0, 1.0)
  grey90       <- rgb(0.9, 0.9, 0.9)
  tv_yellow    <- rgb(1.0, 1.0, 0.2)
  lightmagenta <- rgb(1.0, 0.2, 0.8)
  deepblue     <- rgb(0.25, 0.25, 0.65)
  
  colour_hex   <- c(lightblue, palegreen, yellow, grey60, 
                    grey80, tv_red, tv_green, purpleblue, 
                    grey90, tv_yellow, lightmagenta, deepblue,
                    grey90)
  
  colour_name  <- c("TCRa", "TCRb", "peptide", "MHCa", "MHCb",
                    "CDR1a", "CDR2a", "CDR3a", "FWa",
                    "CDR1b", "CDR2b", "CDR3b", "FWb")
  
  names_in_contact <- unique(as.character(df$Donor_Annotation))
  name_indexes     <- match(sort(names_in_contact), colour_name)
  sub_pallete      <- colour_hex[name_indexes]
  return(sub_pallete)
}

get_start_of_loop <- function(string)
{
  residues <- strsplit(string, '"')
  residues <- residues[[1]]
  residues <- residues[c(3,5,7,9)]
  residues <- gsub("\\[", "", residues)
  residues <- strsplit(residues, ",")
  return(unlist(map(residues, 1)))
}

subtract_CDR <- function(df, list, chain)
{
  #swap the cdrfw and 3 loops arounds
  list     <- list[c(1,2,4,3)]
  loops    <- c("CDR1", "CDR2", "CDR3", "CDR3")
  loops    <- paste(loops, chain, sep = "")
  loops[4] <- paste(loops[4], "fw", sep = "")

  df$Residue <- gsub("A", ".1", df$Residue)
  df$Residue <- gsub("B", ".2", df$Residue)

  list <- gsub("A", ".1", df$Residue)
  list <- gsub("B", ".2", df$Residue)

  for(i in 1:4)
  {
    str   <- loops[i]
    subdf <- df[df$`CDR loop` == str,]
    #print(nrow(subdf))

    if(nrow(df) != 0)
    {
      subdf$Residue <- as.numeric(as.numeric(subdf$Residue) - as.numeric(list[i]))
    }
    else
    {
      next
    }
    if(exists("out") == FALSE)
    {
      out <- subdf
    }
    else
    {
      out <- rbind(out, subdf)
    }
  }
  return(out)
}

get_chain <- function(df, string)
{
  x <- df[grep(string, df$Acceptor_Chain),]
  x <- x[is.na(x$Donor_Annotation) == FALSE,]
  return(x)
}

make_circos <- function(mat, color_vector, color_name, outname, size)
{
  png(outname, 1000, res=120)
  names(color_vector) <- color_name

  circos.par(gap.after = c(rep(5, nrow(mat)-1), 15,
                           rep(5, ncol(mat)-1), 15))
  par(cex = size, mar = c(0, 0, 0, 0))

  chordDiagram(mat, grid.col = color_vector, directional = 2, grid.border = NA,
               annotationTrack = c("grid", "name"))

  dev.off()
  circos.clear()
}

count_contacts <- function(df)
{
  return(df %>%
           group_by(Donor_Annotation, Acceptor_Chain) %>%
           tally)
}

args <- commandArgs(TRUE) #don't need for now

if(length(args) != 4)
{
  print("4 arguments must be supplied")
  print(length(args))
  quit()
}
#read in our files
files <- c(args[1], args[2])

#make outdir
main_dir <- getwd()
sub_dir  <- paste(args[4], "/R_visualisation", sep = "")
dir.create(file.path(main_dir, sub_dir), showWarnings = FALSE)

contacts <- lapply(files, function(i){
  read.table(i, row.names = NULL, sep = "\t")
})

contacts <- do.call(rbind, contacts)

#fix column names being shifted because of rownames being awkward
cols <- colnames(contacts)
cols <- cols[2:length(cols)]

colnames(contacts) <- cols

contacts[,17] <- NULL

#now let's start ploting
lightblue    <- rgb(0.75, 0.75, 1.0)
palegreen    <- rgb(0.65, 0.9, 0.65)
yellow       <- rgb(1.0 , 1.0,	0.0)
grey60       <- rgb(0.6 , 0.6,	0.6)
grey80       <- rgb(0.8 , 0.8,	0.8)
tv_red       <- rgb(1.0, 0.2, 0.2)
tv_green     <- rgb(0.2, 1.0, 0.2)
purpleblue   <- rgb(0.5, 0.0, 1.0)
grey90       <- rgb(0.9, 0.9, 0.9)
tv_yellow    <- rgb(1.0, 1.0, 0.2)
lightmagenta <- rgb(1.0, 0.2, 0.8)
deepblue     <- rgb(0.25, 0.25, 0.65)

colour_hex   <- c(lightblue, palegreen, yellow, grey60,
                  grey80, tv_red, tv_green, purpleblue,
                  grey90, tv_yellow, lightmagenta, deepblue,
                  grey90)
colour_name  <- c("TCRa", "TCRb", "peptide", "MHCa", "MHCb",
                  "CDR1a", "CDR2a", "CDR3a", "FWa",
                  "CDR1b", "CDR2b", "CDR3b", "FWb")



#now let's get our dataframe ready
#first let's get our CDR loops
TCR <- contacts[grep("TCR", contacts$Donor_Chain),]
TCR <- TCR[!is.na(TCR$Donor_Annotation),]

#order as -> 111, 111A, 111B, 112B, 112A, 112
orderings <- c(as.character(seq(1, 110)), c("111", "111A", "111B", "112B", "112A", "112"), as.character(seq(113, 250)))
TCR$Donor_ResNum <- factor(TCR$Donor_ResNum, levels = orderings)

TCR %>%
  group_by(Donor_Annotation, Donor_ResNum, Type) %>%
  tally() %>%
  filter(Donor_Annotation != "") %>%
  ggplot(., aes(x=Donor_ResNum, y=n, fill=Type))+
  geom_bar(stat="identity")+
  facet_grid(Donor_Annotation~., scales = "free")+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  ylab("No. contacts")+
  xlab("residue number")+
  guides(fill=guide_legend(title="Contact force"))+
  theme(plot.title=element_text(hjust = 0.5))

cdr_plot <- paste(sub_dir, "/single_facet.tiff", sep = "")
ggsave(cdr_plot)


TCR %>%
  group_by(Donor_Annotation, Donor_ResNum, Type) %>%
  tally() %>%
  filter(Donor_Annotation != "") %>%
  mutate(chain = replace(Donor_Annotation, 
                         Donor_Annotation == "CDR1a", "alpha")) %>%
  mutate(chain = replace(chain, 
                         chain == "CDR2a", "alpha")) %>%
  mutate(chain = replace(chain, 
                         chain == "CDR3a", "alpha")) %>%
  mutate(chain = replace(chain, 
                         chain == "FWa", "alpha")) %>%
  mutate(chain = replace(chain, 
                         chain == "CDR1b", "beta")) %>%
  mutate(chain = replace(chain, 
                         chain == "CDR2b", "beta")) %>%
  mutate(chain = replace(chain, 
                         chain == "CDR3b", "beta")) %>%
  mutate(chain = replace(chain, 
                         chain == "FWb", "beta")) -> chain_split

chain_split$Donor_Annotation <- gsub('.{1}$', '', chain_split$Donor_Annotation)

ggplot(chain_split, aes(x=Donor_ResNum, y=n, fill=Type))+
  geom_bar(stat="identity")+
  facet_grid(Donor_Annotation~chain, scales = "free")+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  ylab("No. contacts")+
  xlab("residue number")+
  guides(fill=guide_legend(title="Contact force"))+
  theme(legend.title.align = 0.5)

cdr_rev <- paste(sub_dir, "/cdr_loop_force_contact.tiff", sep = "")
ggsave(cdr_rev)

ggplot(chain_split, aes(x=Donor_ResNum, y=n, fill=Type))+
  geom_bar(stat="identity")+
  facet_grid(chain~Donor_Annotation, scales = "free")+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  ylab("No. contacts")+
  xlab("residue number")+
  guides(fill=guide_legend(title="Contact force"))+
  theme(plot.title=element_text(hjust = 0.5))

cdr_facet <- paste(sub_dir, "/cdr_loop_force_contact_rev.tiff", sep = "")
ggsave(cdr_facet)


TCR %>%
  group_by(Donor_Annotation, Type) %>%
  tally() %>%
  filter(Donor_Annotation != "") %>%
  ggplot(., aes(x=Donor_Annotation, y=n, fill=Type))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  ylab("No. contacts")+
  xlab("CDR loop")+
  guides(fill=guide_legend(title="Contact force"))+
  ggtitle("Contact contribution of CDR loop")+
  theme(plot.title=element_text(hjust = 0.5))


force_plot <- paste(sub_dir, "/force_count.tiff", sep = "")
ggsave(force_plot)


###################################################
# pie
##################################################

#get our chains of interest (e.g. cdr -> p and MHC)
to_mhc  <- get_chain(contacts, "MHC")
to_p    <- get_chain(contacts, "peptide")
to_pmhc <- rbind(to_mhc, to_p)

#get our chains of interest (e.g. cdr -> p and MHC)
to_mhc  <- get_chain(contacts, "MHC")
to_p    <- get_chain(contacts, "peptide")
to_pmhc <- rbind(to_mhc, to_p)

to_pmhc %>% group_by(Donor_Annotation, Acceptor_Chain) %>%
  tally() %>%
  filter(Donor_Annotation != "") %>%
  group_by(Acceptor_Chain) %>%
  mutate(freq = sum(n)) %>%
  mutate(group = replace(Acceptor_Chain, Acceptor_Chain == "MHCa", "MHC")) %>% 
  mutate(group = replace(Acceptor_Chain, Acceptor_Chain == "MHCb", "MHC")) %>% 
  mutate(group = replace(Acceptor_Chain, Acceptor_Chain == "peptide", "peptide")) %>%
  ungroup() %>%
  select(-Acceptor_Chain) -> pie_df_split

to_pmhc %>%
  group_by(Donor_Annotation) %>%
  tally() %>%
  mutate(freq = sum(n)) %>%
  filter(Donor_Annotation != "") %>%
  ungroup() -> pie_df_combined

pie_df_combined$group <- "pMHC"

pie_df <- rbind(pie_df_combined, pie_df_split)

#get the palette ready
names_in_contact <- unique(as.character(pie_df$Donor_Annotation))
name_indexes     <- match(sort(names_in_contact), colour_name)
sub_pallete      <- colour_hex[name_indexes]


#non-scaled pie chart
ggplot(pie_df, aes(x=1, y=n, fill=Donor_Annotation, width=1))+
  theme_classic()+
  ggtitle("Contribution of CDR Loop to Contacts")+
  # black border around pie slices
  geom_bar(stat="identity", color='black', position = "fill")+
  # remove black diagonal line from legend
  guides(fill=guide_legend(override.aes=list(colour=NA)))+
  # polar coordinates
  coord_polar(theta='y')+
  # label aesthetics
  theme(axis.ticks=element_blank(),  # the axis ticks
        axis.title=element_blank(),  # the axis labels
        axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels
        axis.text.x=element_blank(),
        plot.title=element_text(hjust = 0.5),
        legend.title.align = 0.5)+
  facet_grid(.~group)+
  scale_fill_manual(values=sub_pallete, name = "CDR Loop")

cdr_pie <- paste(sub_dir, "/cdr_loop_contact_pie.tiff", sep = "")
ggsave(cdr_pie)

print("Making scaled pie")

scaled_pie <- ggplot(pie_df, aes(x=freq/2, y=n, fill=Donor_Annotation, width=freq))+
  theme_classic()+
  ggtitle("Relative Contribution of CDR Loop to Contacts")+
  geom_bar(stat="identity", color='black', position = "fill")+
  guides(fill=guide_legend(override.aes=list(colour=NA)))+
  coord_polar(theta='y')+
  theme(axis.ticks=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        plot.title=element_text(hjust = 0.5),
        legend.title.align = 0.5)+
  facet_grid(.~group)+
  scale_fill_manual(values=sub_pallete, name = "CDR Loop")

cdr_pie_scaled <- paste(sub_dir, "/cdr_loop_contact_scaled_pie.tiff", sep = "")
ggsave(filename = cdr_pie_scaled, plot = scaled_pie)

#now let's have a go at the circos plots
#get rid of ratio
TCR_contacts    <- count_contacts(TCR)
casted_contacts <- dcast(TCR_contacts, Donor_Annotation ~ Acceptor_Chain)

casted_contacts[is.na(casted_contacts)] <- 0

rownames(casted_contacts)        <- casted_contacts$Donor_Annotation
casted_contacts$Donor_Annotation <- NULL
casted_contacts                  <- as.matrix(casted_contacts)
casted_contacts                 <- casted_contacts[rownames(casted_contacts) %in% 
                                                     c("CDR1a","CDR1b","CDR2a", "CDR2b",
                                                       "CDR3a","CDR3b","FWa","FWb"),]


cdr_circos <- paste(sub_dir, "/CDR_circos.png", sep = "")
make_circos(casted_contacts, colour_hex, colour_name, cdr_circos, 0.8)

no_loop_contacts <- TCR %>%
  group_by(Donor_Chain, Acceptor_Chain) %>%
  tally

casted_contacts_2 <- dcast(no_loop_contacts, Donor_Chain ~ Acceptor_Chain)

casted_contacts_2[is.na(casted_contacts_2)] <- 0

rownames(casted_contacts_2)   <- casted_contacts_2$Donor_Chain
casted_contacts_2$Donor_Chain <- NULL
casted_contacts_2             <- as.matrix(casted_contacts_2)

TCR_circos <- paste(sub_dir, "/TCR_circos.png", sep = "")
make_circos(casted_contacts_2, colour_hex, colour_name, TCR_circos, 2.5)
