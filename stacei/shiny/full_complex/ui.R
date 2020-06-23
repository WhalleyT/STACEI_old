pageWithSidebar(
  
  headerPanel("TCR-pMHC Global Contacts"),
  
  sidebarPanel(
    selectInput('facet_row', 'Facet', c(None='.', TCR="Donor_Chain", `CDR Loops`="Donor_Annotation",
                                        pMHC = "Acceptor_Chain", `Interaction Force`="Type")),
    selectInput('colouring', 'Colour', c(None="NULL", TCR="Donor_Chain",`CDR Loops`="Donor_Annotation",
                                         pMHC="Acceptor_Chain", 
                                         `Acceptor Amino Acids`="Acceptor_ResCode", 
                                         `Donor Amino Acids` = "Donor_ResCode",
                                         `Interaction Force`="Type")),
    
    checkboxGroupInput('to_con', 'To', c(TCRa="TCRa", TCRb="TCRb",
                                         `CDR Loops`="CDR Loop")),
    
    checkboxGroupInput('from_con', 'From', c(peptide="peptide", 
                                             MHCa="MHCa", MHCb="MHCb"))
  ),
  
  mainPanel(
    plotOutput('plot')
  )
)

