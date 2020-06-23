setwd("~/Documents/shiny/full_complex")

library("dplyr")
library("tidyr")
library("ggplot2")
library("purrr")
library("RColorBrewer")
library("circlize")
library("reshape2")
library("wesanderson")
sanitize_to <- function(vector)
{
  if(("TCRa" %in% vector || "TCRb" %in% vector) && "CDR Loop" %in% vector)
  {
    return(FALSE)
  }
  else
  {
      return(TRUE)
    }
}

to_column <- function(vector)
{
  if("TCRa" %in% vector || "TCRb" %in% vector)
  {
    return(2)
  }
  else if("CDR Loop" %in% vector)
  {
    return(4)
  }
  else
  {
      return(FALSE)
    }
}

return_subdf <- function(df, to_vec, from_vec)
{
  to_df <- contacts[contacts[,2] %in% to_vec,]
  from_df <- to_df[to_df[,9] %in% from_vec,]
  return(from_df)
}


function(input, output)
{
  dataset <- reactive(
    {
      setwd("/home/tom/Documents/shiny/BSA")
      tcrs <- read.table("TCR_BSA.txt")
      
      colnames(tcrs) <- c("Number", "Chain", "AA", "Residue", "ASA", "BSA", "G", "Gi")
      
      tcrs$Availability <- ((tcrs$ASA - tcrs$BSA) / tcrs$ASA) * 100
      tcrs[["Availability"]][is.na(tcrs[["Availability"]])] <- 0
      
      peptide <- read.table("peptide_BSA.txt")
      peptide$V2 <- as.character(peptide$V2)
      x <- strsplit(peptide$V2, ":")
      peptide_df <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
      
      peptide_df <- cbind(peptide[,1], peptide_df, peptide[,3:7])
      colnames(peptide_df) <- c("Number", "Chain", "AA", "Residue", "ASA", "BSA", "G", "Gi")
      peptide_df$Availability <- ((peptide_df$ASA - peptide_df$BSA) / peptide_df$ASA) * 100
      
      peptide_df[["Availability"]][is.na(peptide_df[["Availability"]])] <- 0
      
      #get chains from fasta
      complexes <- rbind(tcrs, peptide_df)
      
      sequences <- readLines("sequence.fasta")
      for(i in 1:length(sequences))
      {
        line <- sequences[i]
        if(substring(line, 1, 1) ==  ">")
        {
          split <- strsplit(as.character(line), "\\|")
          
          if(split[[1]][3] == "TCRA")
          {
            tcra <- split[[1]][2]
          }
          else if(split[[1]][3] == "TCRB")
          {
            tcrb <- split[[1]][2]
          }
          else if(split[[1]][3] == "peptide")
          {
            peptide <- split[[1]][2]
          }
        }
        
      }
      
      
      complex <- as.character(complexes$Chain)
      complex <- replace(complex, complex == tcra, "TCRa")
      complex <- replace(complex, complex == tcrb, "TCRb")
      complex <- replace(complex, complex == peptide, "peptide")
      complexes$complex <- complex
     
      return(complexes)
      
      #out_cont   <- return_subdf(complexes, input$to_con, input$from_con)

   #   force_count <- out_cont %>%
    #    group_by(Donor_Chain, Donor_ResNum, Acceptor_Chain, Acceptor_ResNum,
    #             Acceptor_ResCode, Donor_ResCode, Type, Donor_Annotation) %>%
       # tally
      #return(force_count)
      
    })
  output$plot <- renderPlot({
  
    
    
    fillvar  <- input$colouring
    facetvar <- input$facet_row
    xaxis    <- input$xaxis
    yaxis    <- input$yaxis
    cols <- colorRampPalette(brewer.pal(9, "Set1"))
    
    
    p <- ggplot(dataset(), aes_string(x="Residue", y=yaxis, fill=fillvar))+
                geom_bar(stat = "identity")+
                xlab("residue")+
                ylab("no. of contacts")+
                theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())+
                scale_fill_manual(values = cols(20))
      
    facets <- paste(facetvar, "~", ".", sep=" ")
    print(facets)
    if(facets != ". ~ .")
    {
      p <- p + facet_grid(facets)
    }
    print(p)
      
  }, height=700)
}



