#!/usr/bin/env python2.7
# This allows us to use the correct display if one is using this script over a ssh session
import matplotlib
matplotlib.use('Agg')

# Our own imports
import bin.new_automated_annotation as annotation
import bin.housekeeping as housekeeping
import bin.containers as classes
import bin.VDJ_usage as vdj
import bin.insertion_aware_VDJ_usage as vdj_out
import bin.renumber_pdb_clean as imgt
import bin.contacts as contacts
import bin.contact_map as con_map
import bin.peptide_pisa as peptide_pisa
import bin.full_pisa as full_pisa
import bin.planar_crossing_angle as crossing_angle
import bin.pymol_cdr_loops as pymol_cdr
import bin.peptide_MHC_visualise as electrostatic
import bin.shape_complementarity as sc

import warnings
import subprocess

from Bio import BiopythonWarning


def main():

    ####################################################################################################################
    """
    This block deals with argument parsing, setting up of paths and variables
    and then the cleaning and annotation of our pdb file.
    """

    print("Collecting parse arguments")
    args, auto = housekeeping.check_parse()
    housekeeping.check_install()


    if args.suppress:
        warnings.simplefilter('ignore', BiopythonWarning)


    print("Assigning classes")
    pdb = classes.PDBStrings(args.infile)


    print("Assigning paths")
    paths = housekeeping.create_paths(pdb.name)


    print("Finding TCR-pMHC chain annotation")
    tcra, tcrb, peptide, mhca, mhcb, mhc_class = annotation.annotate_complex(pdb.file, pdb.filtered, pdb.numbered)
    full_complex = classes.ChainInformation(tcra, tcrb, peptide, mhc_class, mhca, mhcb)


    print("TCRa and TCRb are %s and %s, respectively" % (tcra, tcrb))
    print("MHCa and MHCb are %s and %s respectively" % (mhca, mhcb))
    print("Peptide is %s" % peptide)
    print("MHC is class %i" % mhc_class)

    ####################################################################################################################

    """
    So now we have our basic chain information and a (hopefully) clean PDB file.
    So now we can begin to IMGT number our file
    """
    
    print("Generating FASTA paths")
    anarci_files = classes.AnarciFiles(pdb.name)
    fasta_files = classes.FastaFiles(pdb.name)

    print("Making VDJ assignments")
    vdj.first_annotation_fasta(full_complex, fasta_files.linear, pdb.name, pdb.numbered)

    print("Calling ANARCI")
    vdj.run_anarci(anarci_files.infile, anarci_files.outfile)

    print("Renumbering file to IMGT standards")
    imgt.renumber(anarci_files.outfile, pdb.numbered, full_complex.tcra,
                  full_complex.tcrb, full_complex.mhca, full_complex.mhcb,
                  full_complex.peptide, pdb.imgt)
    
    print("Parsing ANARCI")
    vdj_out.anarci_to_imgt_fasta(pdb.imgt, tcra, tcrb, peptide, mhca, mhcb, anarci_files.outfile, 
                                 pdb.name, fasta_files.linear, fasta_files.annotated)



    ####################################################################################################################

    """
    Complex now should be completely clean! From here we can begin analysis.
    As it is arguably the most important let's start with contacts.
    """
    
    print("Generating contact and sequence paths")
    contact_paths = classes.ContactPaths(pdb.name)
    sequences = classes.LinearSequences(pdb.name)
    tcr_permutations = classes.TCRPermutationContainers()

    # Clean PDB and run both instances of NCONT
    print("Running final pdb cleaning")
    contacts.clean_pdb(pdb.imgt, pdb.name)

    print("Calling NCONT")
    contacts.run_ncont(pdb.name, full_complex.mhca, full_complex.mhcb, full_complex.peptide,
                       full_complex.tcra, full_complex.tcrb, pdb.clean_imgt)

    print("Writing to sequence")
    contacts.pdb_to_sequence(pdb.imgt, full_complex.string, full_complex.mhc_class, pdb.name)

    print("Adding CDR loops")
    contacts.add_cdr_to_sequence(fasta_files.annotated, sequences.seq)

    print("Cleaning and generating TCR to pMHC contacts")
    # TCR -> pMHC contacts
    contacts.clean_contacts(contact_paths.tcr_to_mhc_file, full_complex.string, fasta_files.annotated,
                            args.van_der_waals_distance, args.h_bond_distance, args.s_bond_distance)
    contacts.residue_only_contacts(contact_paths.tcr_to_mhc_clean_file, full_complex.string)
    contacts.annotate_sequence_list(sequences.annotated, contact_paths.tcr_to_mhc_residues)
    
    #compute statistics for CDR loops to major chains
    contacts.stats(contact_paths.tcr_to_mhc_clean_file, pdb.name)


    print("Generating contact maps for TCR to pMHC contacts")

    
    for tcr, pmhc in zip(tcr_permutations.tcr, tcr_permutations.pmhc):
        print(tcr, pmhc)
        con_map.generate_tcr(contact_paths.tcr_to_mhc_list, tcr, pmhc, [], pdb.name)
    for tcr, pmhc, smart in zip(tcr_permutations.tcr_safe, tcr_permutations.pmhc_safe, tcr_permutations.safe_calls):
        con_map.generate_tcr(contact_paths.tcr_to_mhc_list, tcr, pmhc, smart, pdb.name)
    

    print("Cleaning and generating p to MHC contacts")
    # MHC -> peptide contacts
    contacts.clean_contacts(contact_paths.mhc_to_pep_file, full_complex.string, fasta_files.annotated,
                            args.van_der_waals_distance, args.h_bond_distance, args.s_bond_distance)
    contacts.residue_only_contacts(contact_paths.mhc_to_pep_clean_file, full_complex.string)
    contacts.annotate_sequence_list(sequences.annotated, contact_paths.mhc_to_pep_residues)

    print("Generating contact maps for p to MHC")
    
    con_map.generate_mhc(contact_paths.mhc_to_pep_list, full_complex.mhc_class, pdb.name)
    

    ####################################################################################################################

    """
    Now let's call PISA. This will calculate the BSA for the whole TCR-pMHC complex and label it according to CDR loop.
    Then we can do the same, but just for pMHC.
    """
    
    pisa_files = classes.PisaOutputs(pdb.name, full_complex.mhca, full_complex.mhcb, full_complex.peptide,
                                     full_complex.tcra, full_complex.tcrb)

    print("Making a pMHC only PDB file for PISA")
    peptide_pisa.make_pmhc_pdb(pdb.clean_imgt, pdb.pmhc, full_complex.pMHC)

    print("Calling PISA on pMHC complex")
    peptide_pisa.call_pisa(pdb.pmhc, "pMHC_only")
    pepBSA, pepASA = peptide_pisa.extract_pmhc_pisa("pMHC_only", full_complex.peptide, pisa_files.pmhc_chains)

    #Now for the full complex chains
    full_pisa.call_pisa(pdb.imgt, "full_complex")

    for i, j in zip(pisa_files.order, pisa_files.monomers):
        full_pisa.extract_pisa("full_complex", j, full_complex.complex_list, i, full_complex.annotation_dictionary)



    

 ####################################################################################################################

    """
     Pymol based analysis: both visualisation and analysis for crossing angle
    """

    crossing_angle.calculate_and_print(pdb.clean_imgt, fasta_files.annotated, full_complex.mhc_class,
                                       args.ray_trace, full_complex.complex, pdb.name)
    

    
    pymol_cdr.generate(pdb.clean_imgt, fasta_files.annotated, full_complex.mhc_class,
                      full_complex.string, args.ray_trace, pdb.name)
    

    
    
    #buried surface area viz


 ####################################################################################################################

    """
     Check crystal structure validation in pymol; do so with original file so we can 
     keep the CRYST information etc.
    """
    if args.mtz != "None":
        electrostatic.visualise_omit_MHC_only(pdb.file, args.mtz, full_complex.mhc_class,
                                              full_complex.complex, args.ray_trace, pdb.name)

        electrostatic.omit_map(pdb.file, args.mtz, full_complex.mhc_class,
                               full_complex.complex, args.ray_trace, pdb.name)
    else:
        print("Skipping electrostatics")


 ####################################################################################################################
    """
    Run the R code for BSA, and the circos plots/pies (static).
    """
    
    print("Calling R for BSA of peptide")
    subprocess.call("Rscript bin/R/peptide_BSA.R %s" % pisa_files.pmhc_chains, shell=True)

    print("Making Circos plots")

    subprocess.call("Rscript bin/R/circos_and_pie.R %s %s %s %s" % (contact_paths.mhc_to_pep_clean_file,
                                                                    contact_paths.tcr_to_mhc_clean_file,
                                                                    fasta_files.annotated,
                                                                    pdb.name), shell=True)
    
 ###################################################################################################################
    """
    Surface complementarity
    """

    print("Calling SC")
    
    sc.write_SC_pipe(full_complex.mhca, full_complex.mhcb, full_complex.peptide,
                     full_complex.tcra, full_complex.tcrb)

    sc.run_SC(pdb.clean_imgt)
    


 ###################################################################################################################
    """
    clean up and move things to right paths
    """

    print("Done! Now cleaning up")
    housekeeping.clean_namespace(pdb.name, paths, args.infile)


if __name__ == "__main__":
    main()
