import argparse
import sys
import os
import containers
import pymol


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-F', dest='infile', type=str,
                        help='PDB file of a TCR-pMHC complex', required=True)
    parser.add_argument('--mhc', '-Mh', dest='mhc_class', type=int,
                        help='MHC class. if not supplied, then class will be automatically predicted',
                        required=False)
    parser.add_argument('--chains', '-C', dest='chains', type=str,
                        help='Chains of TCR-pMHC complex, in order of TCR alpha, TCR beta,'
                             'peptide, MHC alpha and MHC beta',
                        required=False)
    parser.add_argument('--ray', '-R', dest='ray_trace', action='store_true',
                        help='Flag. If provided, structures will be ray traced in Pymol. This will affect performance.')
    parser.add_argument('--suppress', '-S', dest='suppress', action='store_true',
                        help='Flag. If provided stdout output will be suppressed (inc. CCP4, Pymol and ANARCI)')
    parser.add_argument('--mtz', '-Mt', dest='mtz', required=False, default='ebi', type=str,
                        help='The MTZ file to be analysed, if None is supplied it will skip')
    args = parser.parse_args()
    return args


def check_parse():
    args = _parse_args()

    classes = [1, 2]

    if args.chains is None and args.mhc_class is None:
        return args, True
    elif args.chains is not None and args.mhc_class is not None:
        if args.mhc_class not in classes:
            print "MHC class must be 1 or 2"
            sys.exit()

        if len(args.chains) != 5 and args.chains is not None:
            print "Chains argument must be exactly 5 letters"
            sys.exit()
        return args, False
    else:
        print "Chains argument and MHC argument must be supplied together or not all"
        sys.exit()


def read_file(infile, file_type):
    if infile is None:
        sys.exit('No file to read')
    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
        # print 'Reading file: ' + str(infile)
        return open(infile, "r"), infile.rsplit('.', 1)[0]


def detect_os():
    nixs = ['posix', 'linux']
    if os.name in nixs:
        return "nix"
    else:
        return 'windows'


def create_paths(file_name):
    """creates directories for tidiness (if needed)
       and returns the contacts path for us"""
    if file_name == 'bin':
        print 'Please reconsider calling your PDB file "bin" as it clashes with a directory'
        sys.exit()

    seq_path = file_name + "/sequences"
    contact_path = file_name + "/contacts"
    sc_path = file_name + "/SC"
    pisa_path = file_name + "/pisa"
    xing_path = file_name + "/crossingAngle"
    vis_path = file_name + "/visualisation"
    map_path = file_name + "/maps"
    pdb_path = file_name + "/pdbs"


    if not os.path.exists(file_name):
        # print "Creating Directory " + file_name
        os.makedirs(file_name)

    if not os.path.exists(file_name + "/contacts"):
        # print "Creating Directory " + file_name + "/contacts"
        os.makedirs(contact_path)

    if not os.path.exists(file_name + "/sequences"):
        # print "Creating Directory " + file_name + "/sequences"
        os.makedirs(seq_path)

    if not os.path.exists(file_name + "/pisa"):
        # print "Creating Directory " + file_name + "/pisa"
        os.makedirs(pisa_path)

    if not os.path.exists(file_name + "/SC"):
        # print "Creating Directory " + file_name + "/SC"
        os.makedirs(sc_path)

    if not os.path.exists(xing_path):
        # print "Creating Directory " + file_name + "/crossingAngle"
        os.makedirs(xing_path)

    if not os.path.exists(map_path):
        os.makedirs(map_path)

    if not os.path.exists(vis_path):
        os.makedirs(vis_path)

    if not os.path.exists(pdb_path):
        os.makedirs(pdb_path)

    return containers.Paths(seq_path, contact_path, pisa_path, sc_path, xing_path, map_path, vis_path, pdb_path)


def disable_print(suppress):
    if suppress:
        sys.stdout = open(os.devnull, 'w')


def enable_print():
    sys.stdout = sys.__stdout__


def give_tree(startpath):
    print "Output directory is arranged as follows: \n"
    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * level
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))


def clean_namespace(name):

    # remove uneeded files
    os.remove("clean.fasta")
    os.remove("ANARCI.txt")

    if os.path.exists("ab_contact.txt"):
        os.remove("ab_contact.txt")

    os.remove("peptide_BSA_piped.txt")
    os.remove("Rplots.pdf")
    os.remove("session.txt")

    os.rename(name + "_ANARCI_IMGT_annotated.fasta",
              name + "/" + name + "_ANARCI_IMGT_annotated.fasta")
    os.rename(name + "_MHC_to_pep_contacts.txt",
              name + "/" + name + "_MHC_to_pep_contacts.txt")
    os.rename(name + "_MHC_to_pep_contacts_clean.txt",
              name + "/" + name + "_MHC_to_pep_contacts_clean.txt")
    os.rename(name + "_MHC_to_pep_contacts_clean_residues.txt",
              name + "/" + name + "_MHC_to_pep_contacts_clean_residues.txt")
    os.rename(name + "_MHC_to_pep_contacts_clean_residues_contacts_residues_full.txt",
              name + "/" + name + "_MHC_to_pep_contacts_clean_residues_contacts_residues_full.txt")
    os.rename(name + "_TCR_to_pMHC_contacts.txt",
              name + "/" + name + "_TCR_to_pMHC_contacts.txt")
    os.rename(name + "_TCR_to_pMHC_contacts_clean.txt",
              name + "/" + name + "_TCR_to_pMHC_contacts_clean.txt")
    os.rename(name + "_TCR_to_pMHC_contacts_clean_residues.txt",
              name + "/" + name + "_TCR_to_pMHC_contacts_clean_residues.txt")
    os.rename(name + "_TCR_to_pMHC_contacts_clean_residues_contacts_residues_full.txt",
              name + "/" + name + "_TCR_to_pMHC_contacts_clean_residues_contacts_residues_full.txt")
    os.rename(name + "_clean_numbered_imgt.pdb",
              name + "/" + name + "_clean_numbered_imgt.pdb")
    os.rename(name + "_filtered.pdb", name + "/" + name + "_filtered.pdb")
    os.rename(name + "_mhca_pisa_chains.txt", name + "/" + name + "_mhca_pisa_chains.txt")
    os.rename(name + "_mhcb_pisa_chains.txt", name + "/" + name + "_mhcb_pisa_chains.txt")
    os.rename(name + "_numbered.pdb", name + "/" + name + "_numbered.pdb")
    os.rename(name + "_numbered_imgt.pdb", name + "/" + name + "_numbered_imgt.pdb")
    os.rename(name + "_pMHC_only_pisa_chains.txt",
              name + "/" + name + "_pMHC_only_pisa_chains.txt")
    os.rename(name + "_peptide_pisa_chains.txt",
              name + "/" + name + "_peptide_pisa_chains.txt")
    os.rename(name + "_pmhc.pdb", name + "/" + name + "_pmhc.pdb")
    os.rename(name + "_sequence.txt", name + "/" + name + "_sequence.txt")
    os.rename(name + "_sequence_annot.txt", name + "/" + name + "_sequence_annot.txt")
    os.rename(name + "_statistics.txt", name + "/" + name + "_statistics.txt")
    os.rename(name + "_tcra_pisa_chains.txt", name + "/" + name + "_tcra_pisa_chains.txt")
    os.rename(name + "_tcrb_complex_pisa_chains.txt",
              name + "/" + name + "_tcrb_complex_pisa_chains.txt")
    os.rename(name + "_ANARCI.txt", name + "/" + name + "_ANARCI.txt")
    os.rename("BSA.png", name + "/BSA.png")

    os.rename("sc.txt", name + "/" + "SC/sc.txt")

    os.remove("sc_in.txt")
    os.remove("sc_out.txt")

    # pymol may be left open
    pymol.cmd.quit()
