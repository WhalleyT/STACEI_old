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
    parser.add_argument('--mtz', '-Mt', dest='mtz', required=False, default='None', type=str,
                        help='The MTZ file to be analysed, if None is supplied it will skip')
    args = parser.parse_args()
    return args