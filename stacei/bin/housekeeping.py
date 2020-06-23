import argparse
import sys
import os
import pkg_resources
from . import containers
import glob

from distutils.spawn import find_executable

def _parse_args():
    """
    Parse the command line arguments
    """

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
                        help='Flag. If provided, structures will be ray traced in Pymol. '
                             'This will affect performance at run time but will produce better images.')
    parser.add_argument('--suppress', '-S', dest='suppress', action='store_true',
                        help='Flag. If provided stdout output will be suppressed (inc. CCP4, Pymol and ANARCI)')
    parser.add_argument('--mtz', '-Mt', dest='mtz', required=False, type=str,
                        help='The MTZ file to be analysed, if None (default) is supplied it will skipped',
                        default="None")
    parser.add_argument("--vanderwaals", "-Vdw", dest='van_der_waals_distance', required=False, type=float,
                        default=4.0, help="Distance for an interaction to be considered a contact. If it does not "
                                          "satisfy other criteria it will be marked a van der Waals interaction")
    parser.add_argument("--hydrogenbond", "-Hb", dest="h_bond_distance", required=False, type=float, default=3.5,
                        help="Distance threshold required for a pair of hydrogen bond donors and acceptors to be "
                             "considered to make a hydrogen bond")
    parser.add_argument("--saltbridge", "-Sb", dest='s_bond_distance', required=False, type=float,
                        default=4.0, help="Distance threshold required for a pair of hydrogen bond donors and "
                                          "acceptors to be considered to make a salt bridge")
    args = parser.parse_args()
    return args


def check_parse():
    """
    Sanitize the parse arguments
    """

    args = _parse_args()

    classes = [1, 2]

    if args.chains is None and args.mhc_class is None:
        return args, True
    elif args.chains is not None and args.mhc_class is not None:
        if args.mhc_class not in classes:
            print("MHC class must be 1 or 2")
            sys.exit()

        if len(args.chains) != 5 and args.chains is not None:
            print("Chains argument must be exactly 5 letters")
            sys.exit()
        return args, False
    else:
        print("Chains argument and MHC argument must be supplied together or not all")
        sys.exit()


def read_file(infile, file_type):
    """
    read the file
    """
    if infile is None:
        sys.exit('No file to read')
    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
        # print 'Reading file: ' + str(infile)
        return open(infile, "r"), infile.rsplit('.', 1)[0]


def detect_os():
    """
    Used to deal with faulty calls to
    subprocess, left in just in case
    """
    nixs = ['posix', 'linux']
    if os.name in nixs:
        return "nix"
    else:
        return 'windows'


def create_paths(file_name):
    """creates directories for tidiness (if needed)
       and returns the contacts path for us"""
    if file_name == 'bin':
        print('Please reconsider calling your PDB file "bin" as it clashes with a directory')
        sys.exit()

    seq_path = file_name + "/sequences"
    contact_path = file_name + "/contacts"
    pisa_path = file_name + "/buried_surface"
    sc_path = file_name + "/surface_complementarity"
    xing_path = file_name + "/crossingAngle"
    pdb_path = file_name + "/pdbs"
    vis_path = file_name + "/visualisation"
    session_path = file_name + "/sessions"
    fasta_path = file_name + "/FASTAs"
    elec_path = file_name + "/electrostatics"
    r_path = file_name + "/R_visualisation"

    table_path = contact_path + "/contact_tables"
    mhc_table = table_path + "/MHC_to_pep"
    tcr_table = table_path + "/TCR_to_pMHC"

    paths = [seq_path, contact_path, pisa_path, sc_path, xing_path,
             pdb_path, vis_path, session_path, fasta_path, table_path, mhc_table,
             tcr_table, elec_path, r_path]

    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)

    return containers.Paths(seq_path, contact_path, pisa_path, sc_path, xing_path,
                            pdb_path, vis_path, session_path, fasta_path, elec_path, r_path)


def clean_namespace(name, paths, original):
    """
    Now the script has finished, so time to clean up our files;
    either delete the unwanted ones, or move them to a new home
    """

    # remove uneeded files
    stray_files = ["clean.fasta", "ANARCI.txt", "ab_contact.txt", "session.txt", "peptide_BSA_piped.txt", "Rplots.pdf",
                   "sc_in.txt", "sc_out.txt"]
    for file in stray_files:
        if os.path.exists(file):
            os.remove(file)

    # delete extra fastas
    for file in glob.glob("*.fasta"):
        os.rename(file, paths.fasta_path + "/" + file)

    for file in glob.glob(name + "*.pdb") + glob.glob(paths.crossing_path + "/*.pdb"):
        if file != original:
            os.rename(file, paths.pdb_path + "/" + file.split("/")[-1])

    os.rename("sc.txt", paths.sc_path + "/" + "sc.txt")
    os.rename(name + "_statistics.txt", paths.contact_path + "/contact_tables/" + name + "_statistics.txt")
    os.rename(name + "_ANARCI.txt", paths.sequence_path + "/" + name + "_ANARCI.txt")
    os.rename("BSA.png", paths.r_plots_path + "/BSA.png")


    for file in glob.glob(name + "*pisa_chains*"):
        os.rename(file, paths.pisa_path + "/" + file)

    for file in glob.glob(name + "*MHC_to_pep*"):
        os.rename(file, paths.contact_path + "/contact_tables/MHC_to_pep/" + file)

    for file in glob.glob(name + "*TCR_to_pMHC*"):
        os.rename(file, paths.contact_path + "/contact_tables/TCR_to_pMHC/" + file)

    for file in glob.glob(name + "*sequence*"):
        os.rename(file, paths.sequence_path + "/" + file)

def check_install():
    commands = ["ncont", "pymol", "sc", "ncont", "ANARCI"]

    for comm in commands:
        if find_executable(comm) == False:
            sys.exit("Could not find executable for %s" %comm)


def read_pkg_file(resource_path):
    resource_package = __name__
    print(__name__)
    print(resource_path)
    return pkg_resources.resource_filename(resource_package, resource_path)