#!/usr/bin/env python2.7

# Our own imports
import bin.new_automated_annotation as annotation
import bin.housekeeping as housekeeping
import bin.containers as classes
import bin.VDJ_usage as vdj
import bin.imgt_pdb as imgt
import bin.contacts as contacts
import bin.contact_map as con_map
import bin.new_pisa as pisa
import bin.planar_crossing_angle as crossing_angle
import bin.pymol_cdr_loops as pymol_cdr
import bin.electrostatics as electro

import warnings

from Bio.PDB import PDBExceptions
from Bio import BiopythonWarning


def main():

    ####################################################################################################################
    """
    This block deals with argument parsing, setting up of paths and variables
    and then the cleaning and annotation of our pdb file.
    """

    print "Collecting parse arguments"
    args, auto = housekeeping.check_parse()

    if args.suppress:
        warnings.simplefilter('ignore', PDBExceptions)
        warnings.simplefilter('ignore', BiopythonWarning)

    print "Assigning classes"
    pdb = classes.PDBStrings(args.infile)

    print "Assigning paths"
    paths = housekeeping.create_paths(pdb.name)

    print "Finding TCR-pMHC chain annotation"
    tcra, tcrb, peptide, mhca, mhcb, mhc_class = annotation.annotate_complex(pdb.file, pdb.filtered, pdb.numbered)
    full_complex = classes.ChainInformation(tcra, tcrb, peptide, mhc_class, mhca, mhcb)

    print "TCRa and TCRb are %s and %s, respectively" % (tcra, tcrb)
    print "MHCa and MHCb are %s and %s respectively" % (mhca, mhcb)
    print "Peptide is %s" % peptide
    print "MHC is class %i" % mhc_class

    ####################################################################################################################

    """
    So now we have our basic chain information and a (hopefully) clean PDB file.
    So now we can begin to IMGT number our file
    """

    print "Generating FASTA paths"
    anarci_files = classes.AnarciFiles(pdb.name)
    fasta_files = classes.FastaFiles(pdb.name)

    print "Making VDJ assignments"
    vdj.first_annotation_fasta(full_complex, fasta_files.linear, pdb.name, pdb.numbered)

    print "Calling ANARCI"
    vdj.run_anarci(anarci_files.infile, anarci_files.outfile)

    print "Parsing ANARCI"
    vdj.parse_anarci(full_complex, fasta_files.default, anarci_files.outfile)

    print "Renumbering file to IMGT standards"
    imgt.renumber_ANARCI(fasta_files.annotated, pdb.numbered, full_complex.tcra, full_complex.tcrb,
                         full_complex.peptide, full_complex.mhca, full_complex.mhcb)

    ####################################################################################################################

    """
    Complex now should be completely clean! From here we can begin analysis.
    As it is arguably the most important let's start with contacts.
    """

    print "Generating contact and sequence paths"
    contact_paths = classes.ContactPaths(pdb.name)
    sequences = classes.LinearSequences(pdb.name)
    tcr_permutations = classes.TCRPermutationContainers()

    # Clean PDB and run both instances of NCONT
    print "Running final pdb cleaning"
    contacts.clean_pdb(pdb.imgt, pdb.name)

    print "Calling NCONT"
    contacts.run_ncont(pdb.name, full_complex.mhca, full_complex.mhcb, full_complex.peptide,
                       full_complex.tcra, full_complex.tcrb, pdb.clean_imgt)

    print "Writing to sequence"
    contacts.pdb_to_sequence(pdb.imgt, full_complex.string, full_complex.mhc_class, pdb.name)

    print "Adding CDR loops"
    contacts.add_cdr_to_sequence(fasta_files.annotated, sequences.seq)

    print "Cleaning and generating TCR to pMHC contacts"
    # TCR -> pMHC contacts
    contacts.clean_contacts(contact_paths.tcr_to_mhc_file, full_complex.string, fasta_files.annotated)
    contacts.residue_only_contacts(contact_paths.tcr_to_mhc_clean_file, full_complex.string)
    contacts.annotate_sequence_list(sequences.annotated, contact_paths.tcr_to_mhc_residues)
    contacts.stats(contact_paths.tcr_to_mhc_clean_file, pdb.name)

    """
    print "Generating contact maps for TCR to pMHC contacts"
    for tcr, pmhc in zip(tcr_permutations.tcr, tcr_permutations.pmhc):
        con_map.generate_tcr(contact_paths.tcr_to_mhc_list, tcr, pmhc, [], pdb.name)
    for tcr, pmhc, smart in zip(tcr_permutations.tcr_safe, tcr_permutations.pmhc_safe, tcr_permutations.safe_calls):
        con_map.generate_tcr(contact_paths.tcr_to_mhc_list, tcr, pmhc, smart, pdb.name)
    """

    print "Cleaning and generating p to MHC contacts"
    # MHC -> peptide contacts
    contacts.clean_contacts(contact_paths.mhc_to_pep_file, full_complex.string, fasta_files.annotated)
    contacts.residue_only_contacts(contact_paths.mhc_to_pep_clean_file, full_complex.string)
    contacts.annotate_sequence_list(sequences.annotated, contact_paths.mhc_to_pep_residues)
    contacts.stats(contact_paths.mhc_to_pep_clean_file, pdb.name)

    print "Generating contact maps for p to MHC"
    con_map.generate_mhc(contact_paths.mhc_to_pep_list, full_complex.mhc_class, pdb.name)

    ####################################################################################################################

    """
    Now let's call PISA. This will calculate the BSA for the whole TCR-pMHC complex and label it according to CDR loop.
    Then we can do the same, but just for pMHC.
    """

    pisa_files = classes.PisaOutputs(pdb.name)

    print "Making a pMHC only PDB file for PISA"
    pisa.make_pmhc_pdb(pdb.clean_imgt, pdb.pmhc, full_complex.pMHC)

    print "Calling PISA on pMHC complex"
    pisa.call_pisa(pdb.pmhc, "pMHC_only")
    pepBSA, pepASA = pisa.extract_pisa("pMHC_only", "3", pisa_files.pmhc_chains)

    print "Calling PISA on full complex"
    total_complex_bsa = []


    pisa.call_pisa(pdb.imgt, "full_complex")
    for i,j in zip(pisa_files.order, pisa_files.monomers):
        BSA, ASA = pisa.extract_pisa("full_complex", j, i)

 ####################################################################################################################

    """
     Pymol based analysis: both visualisation and analysis for crossing angle
    """

    """
    crossing_angle.calculate_and_print(pdb.clean_imgt, fasta_files.annotated, full_complex.mhc_class,
                                       args.ray_trace, full_complex.complex)

    pymol_cdr.generate(pdb.clean_imgt, fasta_files.annotated, full_complex.mhc_class,
                      full_complex.string, args.ray_trace, pdb.name)
    """

 ####################################################################################################################

    """
     Check crystal structure validation in pymol
    """
    mtz = None
    if mtz is not None:
        electro.mtz_map()


 ####################################################################################################################

    """
     Check crystal structure validation in pymol
    """




if __name__ == "__main__":
    main()