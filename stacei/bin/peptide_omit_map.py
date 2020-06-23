import os
import pymol
import sys
import argparse
import time
import stacei.bin.data.colourSet as colourSet
import stacei.bin.data.viewSet as viewSet

description = \
    "Generates a peptide omit map by removing the peptide from the PDB model and running a single cycle refmac5 refinement using mtz data either provided via\
    the --MTZ flag or by scraping the EDM entry on PDBe server. The refmac refinement file is then FFT to generate 2FO-1FC and F1=DELFWT PHI=PHDELWT maps. \
    These maps are sent into PyMOL for visualsation of the difference map peaks in the absence of peptide model."


### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--PDB', dest='pdb', type=str, required=True, help='The PDB file to be analysed')
    parser.add_argument('--MTZ', dest='mtz', type=str, required=False, help='The MTZ file to be analysed',
                        default="ebi")
    parser.add_argument('--MHCclass', dest='MHCclass', type=str, required=True,
                        help='The MHC class i.e. class I or II of the complex structure')
    parser.add_argument('--Chains', dest='chains', type=str, required=False,
                        help='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide,TCRa,TCRb', default="ABCDE")
    parser.add_argument('--Ray', dest='ray', type=str, required=False,
                        help='Do you want to render images. Set to False for speed and dev. Default = True',
                        default=True)

    args = parser.parse_args()
    return args


### File loader ###

def readFile(FILE, fileType):
    if FILE == None:
        sys.exit('No file to read')
    if FILE.split('.')[-1].lower() != str(fileType):
        sys.exit("\nFile extension must be of type ." + str(fileType) + "\n")
    else:
        print('Reading file: ' + str(FILE))
        return open(FILE, "r")


def fileToList(inFile):
    allLines = []
    for line in inFile:
        allLines.append(line)
    return allLines


### Pymol tools

def initialisePymol():
    print("\nInitialising pymol...\n")
    pymol.finish_launching(['pymol', '-n'])
    pymol.cmd.reinitialize()
    # set PyMOL parameters
    pymol.cmd.set("ray_shadows", "0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("ray_opaque_background", "0")
    return None


def wait4ray(query):
    print("Waiting for image to render...")
    while not os.path.exists(query):
        time.sleep(1)
    return None


def rayTime(saveas):
    print("Outputting image.. This may take a few seconds..")
    if os.path.exists(saveas):
        print("Removing " + saveas + " as it already exists!")
        os.remove(saveas)
    time.sleep(10)
    pymol.cmd.png(saveas, ray=1, width=3000, height=3000, dpi=300)
    wait4ray(saveas)
    print("Done! " + str(saveas) + " was outputted")


### BODY ###

# Arguments #
def omit_map():
    args = parse_args()
    pdb = args.pdb
    mtz = args.mtz
    MHCclass = args.MHCclass
    chains = args.chains
    ray = args.ray

    pdb_name = pdb.rsplit('.', 1)[0].lower()

    PDBfile = readFile(pdb, "pdb")

    fileName = pdb.rsplit('.', 1)[0]

    # Sort chains
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]


    if mtz != "ebi":
        print("A map.mtz file was provided!", mtz, "will be moved to", fileName + "/maps/" + fileName + ".mtz")
        MTZfile = readFile(mtz, "mtz")
        import shutil

        shutil.copy(mtz, fileName + "/maps/" + fileName + ".mtz")
        MTZfile = readFile(mtz, "mtz")

    if mtz == "ebi":
        if not os.path.exists(fileName + "/maps/" + fileName + ".mtz"):
            print("Downloading map.mtz for entry", pdb_name, "from the PDBe (EBI)")
            import urllib.request, urllib.parse, urllib.error

            urllib.request.urlretrieve("http://www.ebi.ac.uk/pdbe/coordinates/files/" + pdb_name + "_map.mtz",
                               fileName + "/maps/" + fileName + ".mtz")
        else:
            "Did not need to download from ebi as map.mtz already exists"

        if os.path.exists(fileName + "/maps/" + fileName + ".mtz"):
            print("Download successful!")
        else:
            sys.exit(
                "ERROR! Could not scrape the map.mtz file from the EBI server. Ensure input pdb file is named according to pdb entry name or supply loca mtz file")

    mtz = fileName + "/maps/" + fileName + ".mtz"

    ######### OMIT map specific part ##########################

    tempTxt = "rmchain " + peptidechain + "\n" + "END"
    temp = open(fileName + "/pdbs/" + fileName + "_pdbcurPARAM.tmp", "w")
    temp.write(tempTxt)
    temp.close()

    # delete chain
    os.system(
        "pdbcur XYZIN " + pdb + " XYZOUT " + fileName + "/pdbs/" + fileName + "_nopeptide.pdb" + " < " + fileName + "/pdbs/" + fileName + "_pdbcurPARAM.tmp")
    os.remove(fileName + "/pdbs/" + fileName + "_pdbcurPARAM.tmp")

    # Run one cycle of refmac without peptide in the groove
    os.system(
        "refmac5 XYZIN " + fileName + "/pdbs/" + fileName + "_nopeptide.pdb" + " XYZOUT " + fileName + "/pdbs/" + fileName + "_refmac5_omitmap.pdb" + " HKLIN " + fileName + "/maps/" + fileName + ".mtz" + " HKLOUT " + fileName + "/maps/" + fileName + "_refmac5_omitmap.mtz" + " LIBOUT " + fileName + "/maps/" + fileName + "_refmac_omitmap.cif" + " < " + "bin/OMITparam.tmp")

    os.system("mtzdump " + " HKLIN " + fileName + "/maps/" + fileName + "_refmac5_omitmap.mtz END")

    ### This bit is the same as peptideMHCvisualise except the inputs and outputs are different ###

    os.system(
        "fft HKLIN " + fileName + "/maps/" + fileName + "_refmac5_omitmap.mtz" + " MAPOUT " + fileName + "/maps/" + fileName + "_om.map1.tmp" + " < " + "bin/EDMparam1.tmp")
    os.system(
        "mapmask MAPIN " + fileName + "/maps/" + fileName + "_om.map1.tmp" + " MAPOUT " + fileName + "/maps/" + fileName + "_om.map.ccp4" + " XYZIN " + fileName + "/pdbs/" + fileName + "_nopeptide.pdb" + " < " + "bin/EDMparam2.tmp")

    os.system(
        "fft HKLIN " + fileName + "/maps/" + fileName + "_refmac5_omitmap.mtz" + " MAPOUT " + fileName + "/maps/" + fileName + "_om.map3.tmp" + " < " + " bin/EDMparam3.tmp")
    os.system(
        "mapmask MAPIN " + fileName + "/maps/" + fileName + "_om.map3.tmp" + " MAPOUT " + fileName + "/maps/" + fileName + "_om.difference_map.ccp4" + " XYZIN " + fileName + "/pdbs/" + fileName + "_nopeptide.pdb" + " < " + "bin/EDMparam4.tmp")

    os.remove(fileName + "/maps/" + fileName + "_om.map3.tmp")

    ###########################################################


    edmap = fileName + "/maps/" + fileName + "_om.map.ccp4"
    EDMfile = readFile(edmap, "ccp4")

    diffmap = fileName + "/maps/" + fileName + "_om.difference_map.ccp4"
    diffmapFile = readFile(diffmap, "ccp4")

    # Sort chains
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

    # Find the MHC helices
    if MHCclass == "I":
        a1locs = list(range(50, 86))
        MHCa1 = ["MHCa"] + a1locs
        a2locs = list(range(140, 176))
        MHCa2 = ["MHCa"] + a2locs

    if MHCclass == "II":
        a1locs = list(range(46, 78))
        MHCa1 = ["MHCa"] + a1locs
        a2locs = list(range(54, 91))
        MHCa2 = ["MHCb"] + a2locs

    ## Let's get started

    initialisePymol()
    pymol.cmd.reinitialize()

    pymol.cmd.load(pdb, "complex")
    pymol.cmd.load(fileName + "/pdbs/" + fileName + "_refmac5_omitmap.pdb", "omitxyz")
    pymol.cmd.load(edmap, fileName + "_map")
    pymol.cmd.load(diffmap, fileName + "_dmap")

    # align to template
    print("\nAligning file to template...\n")
    pymol.cmd.load("bin/" + MHCclass + "_template.pdb")
    pymol.cmd.align("omitxyz", MHCclass + "_template")
    pymol.cmd.matrix_copy("omitxyz", "complex")
    pymol.cmd.matrix_copy("omitxyz", fileName + "_map")
    pymol.cmd.matrix_copy("omitxyz", fileName + "_dmap")
    pymol.cmd.delete(MHCclass + "_template")
    print("\nAlignment to " + MHCclass + "_template.pdb complete!\n")
    pymol.cmd.delete("omitxyz")

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

    # Select the MHCa helices

    name = "MHCa1"
    locs = '+'.join(str(x) for x in MHCa1[1:])
    pymol.cmd.select("MHCa1", selection=MHCa1[0] + "s and resi " + locs)
    name = "MHCa2"
    locs = '+'.join(str(x) for x in MHCa2[1:])
    pymol.cmd.select("MHCa2", selection=MHCa2[0] + "s and resi " + locs)

    # Make electron density map
    pymol.cmd.map_double(fileName + "_map", -1)
    pymol.cmd.isomesh("p_map_1sigma", fileName + "_map", 1.0, "ps", carve=1.6)
    pymol.cmd.isomesh("p_map_05sigma", fileName + "_map", 0.5, "ps", carve=1.6)
    pymol.cmd.set("mesh_width", 0.5)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.set_view(viewSet.peptideView)

    pymol.cmd.show("mesh", "p_map_1sigma")
    pymol.cmd.color("grey50", "p_map_1sigma")

    ## Photo op here
    pymol.cmd.scene(key="PeptideEdm_omitmap_1sig", action="store")
    if ray == True:
        PeptideEdm1sig = fileName + "/visualisation/" + "omitmap_PeptideEdm1sig.png"
        rayTime(PeptideEdm1sig)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.show("mesh", "p_map_05sigma")
    pymol.cmd.color("grey50", "p_map_05sigma")

    ## Photo op here
    pymol.cmd.scene(key="PeptideEdm_omitmap_05sig", action="store")
    if ray == True:
        PeptideEdm05sig = fileName + "/visualisation/" + "omitmap_PeptideEdm05sig.png"
        rayTime(PeptideEdm05sig)

    # Make difference map

    pymol.cmd.hide("mesh", "all")
    pymol.cmd.isomesh("posdiffmesh", fileName + "_dmap", 3.0, "ps", carve=1.6)
    pymol.cmd.color("green", "posdiffmesh")
    pymol.cmd.show("mesh", "posdiffmesh")
    pymol.cmd.isomesh("negdiffmesh", fileName + "_dmap", -3.0, "ps", carve=1.6)
    pymol.cmd.color("red", "negdiffmesh")
    pymol.cmd.show("mesh", "negdiffmesh")

    ## Photo op here
    pymol.cmd.scene(key="omitmap_differencemap", action="store")
    if ray == True:
        differencemap = fileName + "/visualisation/" + "omitmap_differencemap.png"
        rayTime(differencemap)

    # pMHC helices
    pymol.cmd.hide("mesh", "all")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.color(colourSet.generalcolourSet["MHCa"], "MHCa1")
    pymol.cmd.set("cartoon_transparency", 0.5)

    ## Photo op here
    pymol.cmd.scene(key="MHChelixPeptide1", action="store")
    if ray == True:
        MHChelixPeptide1 = fileName + "/visualisation/" + "omitmap_MHChelixPeptide1.png"
        rayTime(MHChelixPeptide1)

    ## Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    ## Photo op here

    pymol.cmd.scene(key="MHChelixPeptide2", action="store")
    if ray == True:
        MHChelixPeptide2 = fileName + "/visualisation/" + "omitmap_MHChelixPeptide2.png"
        rayTime(MHChelixPeptide2)

    # Save the session
    pymol.cmd.save(fileName + "/visualisation/" + fileName + "_omitmap.pse")

    # Quit pymol
    pymol.cmd.quit()
    print('     ~  End peptideMHCvisualisation_omitmap.py v0.1 BETA  ~')














