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