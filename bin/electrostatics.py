import os
import pymol
import sys
import time

import bin.data.colourSet as colourSet


def read_file(infile, file_type):
    if infile is None:
        sys.exit('No file to read')
    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
        print 'Reading file: ' + str(infile)
        return open(infile, "r")


def file_to_list(infile):
    all_lines = []
    for line in infile:
        all_lines.append(line)
    return all_lines


def init_pymol():
    print "\nInitialising pymol...\n"
    pymol.finish_launching(['pymol', '-qeim'])
    pymol.cmd.reinitialize()

    # set PyMOL parameters
    pymol.cmd.set("ray_shadows", "0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("ray_opaque_background", "0")


def wait4ray(query):
    counter = 0
    while not os.path.exists(query):
        print ("=" * counter) + "| " + str(counter)
        time.sleep(1)
        counter += 1
    return None


def ray_time(saveas):
    print "Outputting image.. This may take a few seconds.."
    if os.path.exists(saveas):
        print "Removing " + saveas + " as it already exists!"
        os.remove(saveas)
    time.sleep(10)
    pymol.cmd.png(saveas, ray=1, width=3000, height=3000, dpi=300)
    wait4ray(saveas)
    print "Done! " + str(saveas) + " was outputted"


def arabic_to_roman(mhc):
    if mhc is 1:
        mhc = "I"
    elif mhc is 2:
        mhc = "II"
    else:
        raise Exception("MHC class is neither I or II")
    return mhc


def mtz_map(pdb, mtz, MHCclass, chains, ray):

    MHCclass = arabic_to_roman(MHCclass)

    # Sort chains
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

    pdb_name = pdb.rsplit('.', 1)[0].lower()

    PDBfile = read_file(pdb, "pdb")

    file_name = pdb.rsplit('.', 1)[0]

    # Make output folder #

    if not os.path.exists(file_name):
        print "Creating Directory " + file_name
        os.makedirs(file_name)

    if not os.path.exists(file_name + "/electrostatics"):
        print "Creating Directory " + file_name + "/electrostatics"
        os.makedirs(file_name + "/electrostatics")

    # Let's get started

    init_pymol()
    pymol.cmd.reinitialize()

    pymol.cmd.load(pdb, "complex")

    # align to template
    print "\nAligning file to template...\n"
    pymol.cmd.load("bin/data/" + MHCclass + "_template.pdb")
    pymol.cmd.align("complex", MHCclass + "_template")
    pymol.cmd.matrix_copy("complex", file_name + "_map")
    pymol.cmd.delete(MHCclass + "_template")
    print "\nAlignment to " + MHCclass + "_template.pdb complete!\n"

    # Make chains objects
    pymol.cmd.select("MHCas", selection="chain " + MHCachain)
    pymol.cmd.select("MHCbs", selection="chain " + MHCbchain)
    pymol.cmd.select("ps", selection="chain " + peptidechain)
    pymol.cmd.select("TCRas", selection="chain " + TCRachain)
    pymol.cmd.select("TCRbs", selection="chain " + TCRbchain)

    pymol.cmd.hide("all")
    pymol.cmd.show("sticks", "ps")
    pymol.cmd.color(colourSet.generalcolourSet["p"], "ps")
    pymol.cmd.util.cnc("ps")