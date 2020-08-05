import os
import glob

class PDBStrings:
    """small class that holds various pdb associated filenames"""

    def __init__(self, arg, outdir):

        #use outdir if specified, else use pdb name
        if outdir != "":

            if outdir.endswith("/") == False:
                outdir = outdir + "/"

            self.name = outdir + arg.rsplit(".")[0].split("/")[-1]
            self.id = arg.rsplit(".")[0].split("/")[-1]
        else:
            self.name = arg.rsplit(".")[0].split("/")[-1]
            self.id = self.name

        self.file = arg
        self.filtered = self.name + "_filtered.pdb"
        self.annotation = self.name + "/" + self.name + "_annotation.txt"
        self.numbered = self.name + "_numbered.pdb"
        self.sc = "sc.pdb"
        self.imgt = self.name + "_numbered_imgt.pdb"
        self.clean_imgt = self.name + "_clean_numbered_imgt.pdb"
        self.pmhc = self.name + "_pmhc.pdb"


class ChainInformation:
    """simple class that holds our information about the
       TCR-pMHC complex"""

    def __init__(self, tcra, tcrb, peptide, mhc_class, mhca, mhcb):
        self.tcra = tcra
        self.tcrb = tcrb
        self.peptide = peptide
        self.mhc_class = mhc_class
        self.mhca = mhca
        self.mhcb = mhcb
        self.complex = mhca + mhcb + peptide + tcra + tcrb
        self.complex_list = list(self.complex)
        if mhc_class == 1:
            self.class_string = "MHC beta chain"
        elif mhc_class == 2:
            self.class_string = "beta-2-microglobulin chain"

        self.pMHC = [peptide, mhca, mhcb]
        self.TCR = [tcra, tcrb]
        self.string = mhca + mhcb + peptide + tcra + tcrb

        self.annotation_dictionary = {tcra: "TCRa",
                                      tcrb: "TCRb",
                                      peptide: "peptide",
                                      mhca: "MHCa",
                                      mhcb: "MHCb"}

    def print_chains(self):
        print("")
        print("PDB file is of class %i" %self.mhc_class)
        print("TCR alpha chain is %s" %self.tcra)
        print("TCR beta chain is %s" %self.tcrb)
        print("MHC alpha chain is %s" %self.mhca)
        print("%s is %s" %(self.class_string, self.mhcb))
        print("Peptide chain is %s" %self.peptide)



class Paths:
    def __init__(self, seq_path, contact_path, pisa_path, sc_path, xing_path,
                 pdb_path, vis_path, session_path, fasta_path, elec_path, r_path, basic_path):
        self.sequence_path = seq_path
        self.contact_path = contact_path
        self.sc_path = sc_path
        self.pisa_path = pisa_path
        self.crossing_path = xing_path
        self.vis_path = vis_path
        self.pdb_path = pdb_path
        self.session_path = session_path
        self.fasta_path = fasta_path
        self.elec_path = elec_path
        self.r_plots_path = r_path
        self.current = os.getcwd()
        self.basic_information = basic_path


class FastaPaths:
    def __init__(self, fasta_path, tcr_ann_fasta_path, cdr_ann_fasta_path, sequence_path, annot_sequence_path):
        self.fasta_path = fasta_path
        self.tcr_ann_fasta_path = tcr_ann_fasta_path
        self.cdr_ann_fasta_path = cdr_ann_fasta_path
        self.sequence_path = sequence_path
        self.annot_sequence_path = annot_sequence_path


class ContactPaths:
    def __init__(self, pdb_name):
        self.mhc_to_pep_file = pdb_name + '_MHC_to_pep_contacts.txt'
        self.tcr_to_mhc_file = pdb_name + '_TCR_to_pMHC_contacts.txt'
        self.mhc_to_pep_clean_file = pdb_name + '_MHC_to_pep_contacts_clean.txt'
        self.tcr_to_mhc_clean_file = pdb_name + '_TCR_to_pMHC_contacts_clean.txt'
        self.tcr_to_mhc_residues = pdb_name + "_TCR_to_pMHC_contacts_clean_residues.txt"
        self.mhc_to_pep_residues = pdb_name + "_MHC_to_pep_contacts_clean_residues.txt"
        self.tcr_to_mhc_list = pdb_name + "_TCR_to_pMHC_contacts_clean_residues_contacts_residues_full.txt"
        self.mhc_to_pep_list =  pdb_name + "_MHC_to_pep_contacts_clean_residues_contacts_residues_full.txt"


class PeptideData:
    def __init__(self, pep_residues, pep_sequence, pepstring, peplen):
        self.residue = pep_residues
        self.sequence = pep_sequence
        self.string = pepstring
        self.length = peplen


class AnarciFiles:
    def __init__(self, outdir, id):
        self.infile =  outdir + "/FASTAs/" + id + ".fasta"
        self.outfile = id + "_ANARCI.txt"


class FastaFiles:
    def __init__(self, outdir, id):
        self.default = outdir + "/FASTAs/" + id +  ".fasta"
        self.linear = outdir + "/FASTAs/" + id + "_linear.fasta"
        self.annotated = id + "_ANARCI_IMGT_annotated.fasta"


class LinearSequences:
    def __init__(self, pdb):
        self.seq = pdb + "_sequence.txt"
        self.annotated = pdb + "_sequence_annot.txt"


class TCRPermutationContainers:
    def __init__(self):
        self.tcr = [["CDR3a"], ["CDR3b"], ["CDR3a", "CDR3b"], ["CDR1a", "CDR2a", "CDR3a"],
                    ["CDR1a", "CDR2a", "FWa", "CDR3a"], ["CDR1b", "CDR2b", "CDR3b"],
                    ["CDR1b", "CDR2b", "FWb", "CDR3b"], ["CDR1b", "CDR2b", "FWb", "CDR3b"],
                    ["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"],
                    ["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"],
                    ["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"]]

        self.pmhc = [["peptide"], ["peptide"], ["peptide"], ["peptide"], ["peptide"],
                     ["peptide"], ["peptide"], ["peptide"], ["peptide"],
                     ["MHCa1", "peptide", "MHCa2"], ["MHCa1", "MHCa2"]]

        self.tcr_safe = [["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"],
                         ["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"],
                         ["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"],
                         ["CDR1b", "CDR2b", "FWb", "CDR3b", "CDR3a", "FWa", "CDR2a", "CDR1a"]]

        self.pmhc_safe = [["peptide"],["MHCa1"], ["MHCa2"], ["MHCa1", "peptide", "MHCa2"]]

        self.safe_calls = [["donors"], ["donors"], ["donors"], ["donors", "acceptors"]]


class PisaOutputs:
    def __init__(self, name, mhca, mhcb, peptide, tcra, tcrb):
        self.pmhc_chains = name + "_pMHC_only_pisa_chains.txt"
        self.complex_chains = name + "_full_complex_pisa_chains.txt"

        self.mhca_chains = name + "_mhca_pisa_chains.txt"
        self.mhcb_chains = name + "_mhcb_pisa_chains.txt"
        self.pept_chains = name + "_peptide_pisa_chains.txt"
        self.tcra_chains = name + "_tcra_pisa_chains.txt"
        self.tcrb_chains = name + "_tcrb_pisa_chains.txt"

        self.order = [self.mhca_chains, self.mhcb_chains, self.pept_chains, self.tcra_chains, self.tcrb_chains]
        self.monomers = [mhca, mhcb, peptide, tcra, tcrb]


class cleanUp:
    def __init__(self):
        print("hi")














