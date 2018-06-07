import pandas as pd

data = pd.read_csv("868_TCR_to_pMHC_contacts_clean.txt", sep="\t")
print data["Donor_Annotation "]
cdrs = data[data["Donor_Annotation "].isin(["CDR3a", "CDR3b"])]

print cdrs