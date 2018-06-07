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
                    "CDR1a", "CDR2a", "CDR3a", "CDR3afw",
                    "CDR1b", "CDR2b", "CDR3b", "CDR3bfw")
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

substract_CDR <- function(df, list, chain)
{
  #swap the cdrfw and 3 loops arounds
  list     <- list[c(1,2,4,3)]
  loops    <- c("CDR1", "CDR2", "CDR3", "CDR3")
  loops    <- paste(loops, chain, sep = "")
  loops[4] <- paste(loops[4], "fw", sep = "")
  
  for(i in 1:4)
  {
    str   <- loops[i]
    subdf <- df[df$`CDR loop` == str,]
    print(nrow(subdf))
    
    if(nrow(df) != 0)
    {
      subdf$Residue <- as.numeric(subdf$Residue - as.numeric(list[i]))
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

make_circos <- function(mat, color_vector, color_name, outname)
{
  
  names(color_vector) <- color_name
  
  circos.par(gap.after = c(rep(5, nrow(mat)-1), 15,
                           rep(5, ncol(mat)-1), 15))
  
  chordDiagram(mat, grid.col = color_vector, directional = 2)
  #ggsave(outname)
  circos.clear()
}

count_contacts <- function(df)
{
  return(df %>%
           group_by(Donor_Annotation, Acceptor_Chain) %>%
           tally)  
}

args <- commandArgs(TRUE) #don't need for now

if(length(args) != 3)
{
  print("3 arguments must be supplied")
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
                  "CDR1a", "CDR2a", "CDR3a", "CDR3afw",
                  "CDR1b", "CDR2b", "CDR3b", "CDR3bfw")



#now let's get our dataframe ready
#first let's get our CDR loops
TCR <- contacts[grep("TCR", contacts$Donor_Chain),]
TCR <- TCR[!is.na(TCR$Donor_Annotation),]

print(head(TCR))

force_count <- TCR %>%
  group_by(Donor_Annotation, Type) %>%
  tally

colnames(force_count) <- c("CDR loop", "Force", "Count")

force_count$Force      <- gsub("HB", "Hydrogen Bond", force_count$Force)
force_count$Force      <- gsub("SB", "Salt Bridge", force_count$Force)
force_count$Force      <- gsub("VW", "Van der Waals", force_count$Force)


ggplot(force_count, aes(x = `CDR loop`, y = Count, fill = Force))+
  geom_bar(stat = "identity")+
  ggtitle("Contact Contribution of CDR loops")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title.align = 0.5)+
  scale_fill_brewer(palette = "Set2")
force_count_out <- paste(sub_dir, "/force_count.tiff", sep = "")
ggsave(force_count_out)

force_count_residues <- TCR %>%
  group_by(Donor_Annotation, Donor_ResNum, Type) %>%
  tally

colnames(force_count_residues) <- c("CDR loop", "Residue", "Force", "Count")


#get start residues of fasta
annotated_fasta <- readLines(args[3])
tcra_line       <- annotated_fasta[7]
tcrb_line       <- annotated_fasta[9]

tcra_cdr_starts <- get_start_of_loop(tcra_line)
tcrb_cdr_starts <- get_start_of_loop(tcrb_line)


x <- substract_CDR(as.data.frame(force_count_residues), tcra_cdr_starts, "a")
y <- substract_CDR(as.data.frame(force_count_residues), tcrb_cdr_starts, "b")
z <- rbind(x,y)

colnames(z) <- c("CDR Loop", "Residue", "Force", "Count")

ggplot(z, aes(x = Residue, y = Count, fill = Force))+
  geom_bar(stat = "identity")+
  facet_grid(`CDR Loop` ~ .)+
  scale_x_continuous(breaks = seq(0, max(z$Count), 1))

single_facet <- paste(sub_dir, "/single_facet.tiff", sep = "")
ggsave(single_facet)

x <- z[z$`CDR Loop` %in% c("CDR1a", "CDR2a", "CDR3a", "CDR3afw"),]
y <- z[z$`CDR Loop` %in% c("CDR1b", "CDR2b", "CDR3b", "CDR3bfw"),]
x$Chain <- "alpha"
y$Chain <- "beta"

zz <- rbind(x,y)

zz$`CDR Loop` <- gsub('a', '', zz$`CDR Loop`)
zz$`CDR Loop` <- gsub('b', '', zz$`CDR Loop`)
zz$`CDR Loop` <- gsub('fw', '', zz$`CDR Loop`)
zz$Force      <- gsub("HB", "Hydrogen Bond", zz$Force)
zz$Force      <- gsub("SB", "Salt Bridge", zz$Force)
zz$Force      <- gsub("VW", "Van der Waals", zz$Force)

ggplot(zz, aes(x = Residue, y = Count, fill = Force))+
  geom_bar(stat = "identity")+
  facet_grid(`CDR Loop` ~ Chain)+
  scale_x_continuous(breaks = seq(0, max(zz$Count), 1))+
  ggtitle("Contact Contribution of CDR loops")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title.align = 0.5)+
  scale_fill_brewer(palette = "Set2")

force_count_double <- paste(sub_dir, "/cdr_loop_force_contact.tiff", sep = "")
ggsave(force_count_double)



ggplot(zz, aes(x = Residue, y = Count, fill = Force))+
  geom_bar(stat = "identity")+
  facet_grid(Chain ~ `CDR Loop`)+
  scale_x_continuous(breaks = seq(0, max(zz$Count), 1))

force_count_double_rev <- paste(sub_dir, "/cdr_loop_force_contact_rev.tiff", sep = "")
ggsave(force_count_double_rev)

################
#now pie charts#
################

#get our chains of interest (e.g. cdr -> p and MHC)
to_mhc  <- get_chain(contacts, "MHC")
to_p    <- get_chain(contacts, "peptide")
to_pmhc <- rbind(to_mhc, to_p)

#we can keep above if we want to compare different chains etc.
to_p_count    <- plyr::count(to_p, "Donor_Annotation")
to_mhc_count  <- plyr::count(to_mhc, "Donor_Annotation")
to_pmhc_count <- plyr::count(to_pmhc, "Donor_Annotation")

to_pmhc_count$Source <- "pMHC"
to_mhc_count$Source  <- "MHC"
to_p_count$Source    <- "peptide"

all_mhc_counts <- rbind(to_p_count, to_mhc_count, to_pmhc_count)

total_size <- sum(all_mhc_counts$freq)

to_pmhc_count$Size   <- sum(to_pmhc_count$freq) / total_size
to_mhc_count$Size    <- sum(to_mhc_count$freq) / total_size
to_p_count$Size      <- sum(to_p_count$freq) / total_size

all_mhc_counts <- rbind(to_p_count, to_mhc_count, to_pmhc_count)

names_in_contact <- unique(as.character(all_mhc_counts$Donor_Annotation))
name_indexes     <- match(sort(names_in_contact), colour_name)
sub_pallete      <- colour_hex[name_indexes]

#so there's no way that seems to properly let me fit the pallete
#to a given string. let's try and order our pallete based on what's
#in the mhc dataframe

names_in_contact <- unique(as.character(all_mhc_counts$Donor_Annotation))
name_indexes     <- match(sort(names_in_contact), colour_name)
sub_pallete      <- colour_hex[name_indexes]

#now let's make value for relative size of pie


ggplot(all_mhc_counts, aes(x=Size/2, y=freq, fill=Donor_Annotation,
                           width=Size))+
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
  facet_grid(.~Source)+
  scale_fill_manual(values=sub_pallete, name = "CDR Loop")

cdr_pie <- paste(sub_dir, "/cdr_loop_contact_pie.tiff", sep = "")
ggsave(cdr_pie)


#now let's have a go at the circos plots
#get rid of ratio
TCR_contacts    <- count_contacts(TCR)
casted_contacts <- dcast(TCR_contacts, Donor_Annotation ~ Acceptor_Chain)

casted_contacts[is.na(casted_contacts)] <- 0

rownames(casted_contacts)        <- casted_contacts$Donor_Annotation
casted_contacts$Donor_Annotation <- NULL
casted_contacts                  <- as.matrix(casted_contacts)



cdr_circos <- paste(sub_dir, "/cdr_loop_circos.tiff", sep = "")
make_circos(casted_contacts, colour_hex, colour_name, cdr_circos)

no_loop_contacts <- TCR %>%
  group_by(Donor_Chain, Acceptor_Chain) %>%
  tally

casted_contacts_2 <- dcast(no_loop_contacts, Donor_Chain ~ Acceptor_Chain)

casted_contacts_2[is.na(casted_contacts_2)] <- 0

rownames(casted_contacts_2)   <- casted_contacts_2$Donor_Chain
casted_contacts_2$Donor_Chain <- NULL
casted_contacts_2             <- as.matrix(casted_contacts_2)

TCR_circos <- paste(sub_dir, "/TCR_circos.tiff", sep = "")
make_circos(casted_contacts_2, colour_hex, colour_name, TCR_circos)
