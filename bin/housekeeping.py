import argparse
import sys
import os
import containers


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-F', dest='infile', type=str,
                        help='PDB file of a TCR-pMHC complex', required=True)
    parser.add_argument('--mhc', '-M', dest='mhc_class', type=int,
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

    return containers.Paths(seq_path, contact_path, pisa_path, sc_path, xing_path)


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
