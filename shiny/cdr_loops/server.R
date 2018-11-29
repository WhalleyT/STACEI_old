setwd("~/Documents/shiny/cdr_loops")

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
      tcr <- read.table("tcr.txt", row.names = NULL, sep = "\t")
      mhc <- read.table("mhc.txt", row.names = NULL, sep = "\t")
      contacts <- rbind(tcr,mhc)
      cols <- colnames(contacts)
      cols <- cols[2:length(cols)]
      colnames(contacts) <- cols
      
      correct_to <- sanitize_to(input$to_con)
      to_col     <- to_column(input$to_con)
      out_cont   <- return_subdf(contacts, input$to_con, input$from_con)
      
      #Get only CDR related columns
      CDR_terms <- c("CDR1a", "CDR2a", "CDR3a", "FWa",
                     "CDR1b", "CDR2b", "CDR3b", "FWb")
      
      contacts <- out_cont[out_cont$Donor_Annotation %in% CDR_terms,]
      
      contact_groups <- contacts %>%
        group_by(Donor_Chain, Donor_ResNum, Acceptor_Chain, Acceptor_ResNum,
                 Acceptor_ResCode, Donor_ResCode, Type, Donor_Annotation) %>%
        tally
      return(contact_groups)
      
    })
  output$plot <- renderPlot({
    
    
    
    fillvar  <- input$colouring
    facetvar <- input$facet_row
    cols <- colorRampPalette(brewer.pal(9, "Set1"))
    
    print(dataset)
    
    p <- ggplot(dataset(), aes_string(x="Donor_Annotation", y="n", fill=fillvar))+
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



