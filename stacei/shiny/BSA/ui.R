pageWithSidebar(
  
  headerPanel("Buried Surface Area"),
  
  sidebarPanel(
    selectInput('facet_row', 'Facet', c(None='.', Complex="complex")),
    selectInput('colouring', 'Colour', c(None="NULL", 
                                         `Complex`="complex", 
                                         `Amino Acids` = "AA")),
    
    checkboxGroupInput('to_con', 'To', c(TCRa="TCRa", TCRb="TCRb",
                                         `Peptide`="peptide")),
    selectizeInput('yaxis', 'Y axis', c(Availability="Availability",
                                        `Buried Surface Area` ="BSA",
                                        `Available Surface Area` = "ASA"))
    
    #checkboxGroupInput('from_con', 'From', c(peptide="peptide", 
                                             #MHCa="MHCa", MHCb="MHCb")),
    #selectizeInput('xaxis', 'X axis', c(Residue="Residue", "")
  ),
  
  mainPanel(
    plotOutput('plot')
  )
)

