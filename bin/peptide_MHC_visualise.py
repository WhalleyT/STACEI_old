import os
import pymol
import sys
import time
import shutil
import urllib.request, urllib.parse, urllib.error

import bin.data.colourSet as colourSet
import bin.data.viewSet as viewSet


def empty_file(filename):
    if filename is None:
        sys.exit('No file to read')


def read_file(filename, file_type):
    empty_file(filename)

    if filename.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
        print('Reading file: ' + str(filename))
        return open(filename, "r")


def file_to_list(infile):
    allLines = []
    for line in infile:
        allLines.append(line)
    return allLines


def initialisePymol():
    print("\nInitialising pymol...\n")
    pymol.finish_launching(['pymol', '-nqc'])
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


def ray_tracer(saveas, tracing):
    print("Outputting image.. This may take a few seconds..")
    if os.path.exists(saveas):
        print("Removing " + saveas + " as it already exists!")
        os.remove(saveas)
    time.sleep(10)
    pymol.cmd.png(saveas, ray=tracing, width=3000, height=3000, dpi=300)
    wait4ray(saveas)
    print("Done! " + str(saveas) + " was outputted")


def visualise_omit_MHC_only(pdb, mtz, MHCclass, chains, ray, file_name):
    #directly fetch our pdb file name in case "EBI" has been specified and outdir as outdir flag with override
    pdb_name = pdb.rsplit(".")[0].split("/")[-1].lower()

    prefix = file_name + "/electrostatics/" + pdb_name
    if MHCclass == 1:
        MHCclass = "I"
    elif MHCclass == 2:
        MHCclass = "II"


    # Make output folder #
    print(file_name)

    if not os.path.exists(file_name):
        print("Creating Directory " + file_name)
        os.makedirs(file_name)

    if not os.path.exists(file_name + "/visualisation"):
        print("Creating Directory " + file_name + "/visualisation")
        os.makedirs(file_name + "/visualisation")

    if not os.path.exists(file_name + "/maps"):
        print("Creating Directory " + file_name + "/maps")
        print("Creating Directory " + file_name + "/maps")
        os.makedirs(file_name + "/maps")

    if mtz != "ebi":
        print("A map.mtz file was provided!", mtz, "will be moved to", prefix + ".mtz")

        shutil.copy(mtz, prefix + ".mtz")

    if mtz == "ebi":
        if not os.path.exists(prefix + ".mtz"):
            print("Downloading map.mtz for entry", pdb_name, "from the PDBe (EBI)")

            try:
                url = "http://www.ebi.ac.uk/pdbe/coordinates/files/" + pdb_name + "_map.mtz"
                destination = prefix + ".mtz"
                print(url, destination)
                urllib.request.urlretrieve(url, destination)
            except:
                print("Could not retrieve url. Please try again, making sure you either supply a file," \
                      " or make sure your file shares its name with one on PDB")
                # quit early, get rid of pymol
                pymol.cmd.quit()
                sys.exit()


        else:
            "Did not need to download from ebi as map.mtz already exists"

        if os.path.exists(prefix + ".mtz"):
            print("Download successful!")

    mtz = prefix + ".mtz"

    # Use CCP4 to generate map

    os.system("fft HKLIN " + mtz + " MAPOUT " + file_name + "/electrostatics/" +
              pdb_name + ".map1.tmp" + " < " + "bin/data/EDMparam1.tmp")
    os.system("mapmask MAPIN " + prefix + ".map1.tmp" +
              " MAPOUT " + prefix + ".map.ccp4" + " XYZIN " + pdb + " < " + "bin/data/EDMparam2.tmp")
    os.system("fft HKLIN " + mtz + " MAPOUT " + file_name + "/electrostatics/" +
              pdb_name + ".map3.tmp" + " < " + " bin/data/EDMparam3.tmp")
    os.system("mapmask MAPIN " + prefix + ".map3.tmp" +
              " MAPOUT " + prefix + ".difference_map.ccp4" +
              " XYZIN " + pdb + " < " + "bin/data/EDMparam4.tmp")

    edmap = prefix + ".map.ccp4"

    diffmap = prefix + ".difference_map.ccp4"

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

    # Let's get started

    initialisePymol()
    pymol.cmd.reinitialize()

    pymol.cmd.load(pdb, "complex")
    pymol.cmd.load(edmap, pdb_name + "_map")
    pymol.cmd.load(diffmap, pdb_name + "_dmap")

    # align to template
    print("\nAligning file to template...\n")
    pymol.cmd.load("bin/data/" + MHCclass + "_template.pdb")
    pymol.cmd.align("complex", MHCclass + "_template")
    pymol.cmd.matrix_copy("complex", pdb_name + "_map")
    pymol.cmd.delete(MHCclass + "_template")
    print("\nAlignment to " + MHCclass + "_template.pdb  complete!\n")

    # Make chains objects
    pymol.cmd.select("MHCas", selection="chain " + MHCachain)
    pymol.cmd.select("MHCbs", selection="chain " + MHCbchain)
    pymol.cmd.select("ps", selection="chain " + peptidechain)
    pymol.cmd.select("TCRas", selection="chain " + TCRachain)
    pymol.cmd.select("TCRbs", selection="chain " + TCRbchain)

    pymol.cmd.hide("all")
    pymol.cmd.show("sticks", "ps")
    pymol.cmd.color(colourSet.generalColourSet["p"], "ps")
    pymol.cmd.util.cnc("ps")

    # Select the MHCa helices
    #MHCa1, MHCa2 = None, None

    locs = '+'.join(str(x) for x in MHCa1[1:])
    pymol.cmd.select("MHCa1", selection=MHCa1[0] + "s and resi " + locs)
    locs = '+'.join(str(x) for x in MHCa2[1:])
    pymol.cmd.select("MHCa2", selection=MHCa2[0] + "s and resi " + locs)

    # Make electron density map
    pymol.cmd.map_double(file_name + "_map", -1)
    pymol.cmd.isomesh("p_map_1sigma", pdb_name + "_map", 1.0, "ps", carve=1.6)
    pymol.cmd.isomesh("p_map_05sigma", pdb_name + "_map", 0.5, "ps", carve=1.6)
    pymol.cmd.set("mesh_width", 0.5)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.set_view(viewSet.peptideView)

    pymol.cmd.show("mesh", "p_map_1sigma")
    pymol.cmd.color("grey50", "p_map_1sigma")

    # Photo op here
    pymol.cmd.scene(key="PeptideEdm1sig", action="store")
    PeptideEdm1sig = file_name + "/visualisation/" + pdb_name + "_PeptideEdm1sig.png"

    if ray:
        ray_tracer(PeptideEdm1sig, 1)
    else:
        ray_tracer(PeptideEdm1sig, 0)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.show("mesh", "p_map_05sigma")
    pymol.cmd.color("grey50", "p_map_05sigma")

    # Photo op here
    pymol.cmd.scene(key="PeptideEdm05sig", action="store")
    PeptideEdm05sig = file_name + "/visualisation/" + pdb_name + "_PeptideEdm05sig.png"

    if ray:
        ray_tracer(PeptideEdm05sig, 1)
    else:
        ray_tracer(PeptideEdm05sig, 0)

    # Make difference map

    pymol.cmd.hide("mesh", "all")
    pymol.cmd.isomesh("posdiffmesh", pdb_name + "_dmap", 3.0, "ps", carve=1.6)
    pymol.cmd.color("green", "posdiffmesh")
    pymol.cmd.show("mesh", "posdiffmesh")
    pymol.cmd.isomesh("negdiffmesh", pdb_name + "_dmap", -3.0, "ps", carve=1.6)
    pymol.cmd.color("red", "negdiffmesh")
    pymol.cmd.show("mesh", "negdiffmesh")

    # Photo op here
    pymol.cmd.scene(key="differencemap", action="store")
    differencemap = file_name + "/visualisation/" + pdb_name + "_differencemap.png"

    if ray:
        ray_tracer(differencemap, 1)
    else:
        ray_tracer(differencemap, 0)

    # pMHC helices
    pymol.cmd.hide("mesh", "all")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.color(colourSet.generalColourSet["MHCa"], "MHCa1")
    pymol.cmd.set("cartoon_transparency", 0.5)

    # Photo op here
    pymol.cmd.scene(key="MHChelixPeptide1", action="store")
    MHChelixPeptide1 = file_name + "/visualisation/" + pdb_name + "_MHChelixPeptid1e.png"

    if ray:
        ray_tracer(MHChelixPeptide1, 1)
    else:
        ray_tracer(MHChelixPeptide1, 0)

    # Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    # Photo op here

    pymol.cmd.scene(key="MHChelixPeptide2", action="store")
    MHChelixPeptide2 = file_name + "/visualisation/" + pdb_name + "_MHChelixPeptide2.png"

    if ray:
        ray_tracer(MHChelixPeptide2, 1)
    else:
        ray_tracer(MHChelixPeptide2, 0)

    # Save the session
    pymol.cmd.save(file_name + "/sessions/" + pdb_name + "_peptideMHCvis.pse")

    # Quit pymol
    #pymol.cmd.quit()
    print('     ~  End peptideMHCvisualisation.py v0.1 BETA  ~')


def omit_map(pdb, mtz, MHCclass, chains, ray, file_name):
    pdb_name = pdb.rsplit(".")[0].split("/")[-1].lower()
    prefix = file_name + "/electrostatics/" + pdb_name

    if MHCclass == 1:
        MHCclass = "I"
    elif MHCclass == 2:
        MHCclass = "II"

    # Sort chains
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

    # Make output folder #

    if not os.path.exists(file_name):
        print("Creating Directory " + file_name)
        os.makedirs(file_name)

    if not os.path.exists(file_name + "/visualisation"):
        print("Creating Directory " + file_name + "/visualisation")
        os.makedirs(file_name + "/visualisation")

    if not os.path.exists(file_name + "/maps"):
        print("Creating Directory " + file_name + "/maps")
        os.makedirs(file_name + "/maps")

    if not os.path.exists(file_name + "/pdbs"):
        print("Creating Directory " + file_name + "/pdbs")
        os.makedirs(file_name + "/pdbs")

    if mtz != "ebi":
        print("A map.mtz file was provided!", mtz, "will be moved to", prefix + ".mtz")

        shutil.copy(mtz, prefix + ".mtz")

    if mtz == "ebi":
        if not os.path.exists(prefix + ".mtz"):
            print("Downloading map.mtz for entry", pdb_name, "from the PDBe (EBI)")

            try:
                urllib.request.urlretrieve("http://www.ebi.ac.uk/pdbe/coordinates/files/" + pdb_name + "_map.mtz",
                               prefix + ".mtz")
            except:
                print("Could not retrieve url. Please try again, making sure you either supply a file," \
                      " or your file shares its name with one on PDB")
                # quit early, get rid of pymol
                pymol.cmd.quit()
                sys.exit()


    tempTxt = "delchain " + peptidechain + "\n" + "END"
    temp = open(file_name + "/pdbs/" + pdb_name + "_pdbcurPARAM.tmp", "w")
    temp.write(tempTxt)
    temp.close()

    # delete chain
    print("Deleting chain for %s" % pdb)
    os.system("pdbcur XYZIN " + pdb + " XYZOUT " + file_name + "/pdbs/" + pdb_name +
              "_nopeptide.pdb" + " < " + file_name + "/pdbs/" + pdb_name + "_pdbcurPARAM.tmp")

    # Run one cycle of refmac without peptide in the groove
    print("Running refmac with peptideless complex")
    os.system("refmac5 XYZIN " + file_name + "/pdbs/" + pdb_name + "_nopeptide.pdb" + " XYZOUT " +
              file_name + "/pdbs/" + pdb_name + "_refmac5_omitmap.pdb" + " HKLIN " + file_name + "/electrostatics/"
              + pdb_name + ".mtz" + " HKLOUT " + prefix + "_refmac5_omitmap.mtz" + " LIBOUT "
              + prefix + "_refmac_omitmap.cif" + " < " + "bin/data/OMITparam.tmp")


    os.system("mtzdump  HKLIN " + prefix + "_refmac5_omitmap.mtz END")


    # This bit is the same as peptideMHCvisualise except the inputs and outputs are different
    print("Running fast fourier transform with parameters 1")
    os.system("fft HKLIN " + prefix + "_refmac5_omitmap.mtz" + " MAPOUT " + file_name +
              "/electrostatics/" + pdb_name + "_om.map1.tmp" + " < " + "bin/data/EDMparam1.tmp")
    print("Running mapmask with parameters 1")
    os.system("mapmask MAPIN " + prefix + "_om.map1.tmp" +
              " MAPOUT " + prefix + "_om.map.ccp4" + " XYZIN "
              + file_name + "/pdbs/" + pdb_name + "_nopeptide.pdb" + " < " + "bin/data/EDMparam2.tmp")
    print("Running fast fourier transform with parameters 2")
    os.system("fft HKLIN " + prefix + "_refmac5_omitmap.mtz"
              + " MAPOUT " + prefix + "_om.map3.tmp" + " < " + " bin/data/EDMparam3.tmp")
    print("Running mapmask with parameters 2")
    os.system("mapmask MAPIN " + prefix + "_om.map3.tmp" +
              " MAPOUT " + prefix + "_om.difference_map.ccp4"
              + " XYZIN " + file_name + "/pdbs/" + pdb_name + "_nopeptide.pdb" + " < " + "bin/data/EDMparam4.tmp")

    os.remove(prefix + "_om.map3.tmp")
    os.remove(prefix + ".map1.tmp")
    os.remove(prefix + ".map3.tmp")
    os.remove(file_name + "/pdbs/" + pdb_name + "_pdbcurPARAM.tmp")

    edmap = prefix + "_om.map.ccp4"

    diffmap = prefix + "_om.difference_map.ccp4"

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

    # Let's get started

    initialisePymol()
    pymol.cmd.reinitialize()

    pymol.cmd.load(pdb, "complex")
    pymol.cmd.load(file_name + "/pdbs/" + pdb_name + "_refmac5_omitmap.pdb", "omitxyz")
    pymol.cmd.load(edmap, pdb_name + "_map")
    pymol.cmd.load(diffmap, pdb_name + "_dmap")

    # align to template
    print("\nAligning file to template...\n")
    pymol.cmd.load("bin/data/" + MHCclass + "_template.pdb")
    pymol.cmd.align("omitxyz", MHCclass + "_template")
    pymol.cmd.matrix_copy("omitxyz", "complex")
    pymol.cmd.matrix_copy("omitxyz", pdb_name + "_map")
    pymol.cmd.matrix_copy("omitxyz", pdb_name + "_dmap")
    pymol.cmd.delete(MHCclass + "_template")
    print("\nAlignment to " + MHCclass + "_template.pdb  complete!\n")
    pymol.cmd.delete("omitxyz")

    # Make chains objects
    pymol.cmd.select("MHCas", selection="chain " + MHCachain)
    pymol.cmd.select("MHCbs", selection="chain " + MHCbchain)
    pymol.cmd.select("ps", selection="chain " + peptidechain)
    pymol.cmd.select("TCRas", selection="chain " + TCRachain)
    pymol.cmd.select("TCRbs", selection="chain " + TCRbchain)

    pymol.cmd.hide("all")
    pymol.cmd.show("sticks", "ps")
    pymol.cmd.color(colourSet.generalColourSet["p"], "ps")
    pymol.cmd.util.cnc("ps")

    # Select the MHCa helices
    #MHCa1, MHCa2 = None, None

    locs = '+'.join(str(x) for x in MHCa1[1:])
    pymol.cmd.select("MHCa1", selection=MHCa1[0] + "s and resi " + locs)
    locs = '+'.join(str(x) for x in MHCa2[1:])
    pymol.cmd.select("MHCa2", selection=MHCa2[0] + "s and resi " + locs)

    # Make electron density map
    pymol.cmd.map_double(file_name + "_map", -1)
    pymol.cmd.isomesh("p_map_1sigma", pdb_name + "_map", 1.0, "ps", carve=1.6)
    pymol.cmd.isomesh("p_map_05sigma", pdb_name + "_map", 0.5, "ps", carve=1.6)
    pymol.cmd.set("mesh_width", 0.5)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.set_view(viewSet.peptideView)

    pymol.cmd.show("mesh", "p_map_1sigma")
    pymol.cmd.color("grey50", "p_map_1sigma")

    # Photo op here
    pymol.cmd.scene(key="PeptideEdm_omitmap_1sig", action="store")
    PeptideEdm1sig = file_name + "/visualisation/" + pdb_name + "_omitmap_PeptideEdm1sig.png"

    if ray:
        ray_tracer(PeptideEdm1sig, 1)
    else:
        ray_tracer(PeptideEdm1sig, 0)

    pymol.cmd.hide("mesh", "all")

    pymol.cmd.show("mesh", "p_map_05sigma")
    pymol.cmd.color("grey50", "p_map_05sigma")

    # Photo op here
    pymol.cmd.scene(key="PeptideEdm_omitmap_05sig", action="store")
    PeptideEdm05sig = file_name + "/visualisation/" + pdb_name + "_omitmap_PeptideEdm05sig.png"

    if ray:
        ray_tracer(PeptideEdm05sig, 1)
    else:
        ray_tracer(PeptideEdm05sig, 0)


    pymol.cmd.hide("mesh", "all")
    pymol.cmd.isomesh("posdiffmesh", pdb_name + "_dmap", 3.0, "ps", carve=1.6)
    pymol.cmd.color("green", "posdiffmesh")
    pymol.cmd.show("mesh", "posdiffmesh")
    pymol.cmd.isomesh("negdiffmesh", pdb_name + "_dmap", -3.0, "ps", carve=1.6)
    pymol.cmd.color("red", "negdiffmesh")
    pymol.cmd.show("mesh", "negdiffmesh")

    # Photo op here
    pymol.cmd.scene(key="omitmap_differencemap", action="store")
    differencemap = file_name + "/visualisation/" + pdb_name + "+omitmap_differencemap.png"

    if ray:
        ray_tracer(differencemap, 1)
    else:
        ray_tracer(differencemap, 0)

    # pMHC helices
    pymol.cmd.hide("mesh", "all")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.color(colourSet.generalColourSet["MHCa"], "MHCa1")
    pymol.cmd.set("cartoon_transparency", 0.5)

    # Photo op here
    pymol.cmd.scene(key="MHChelixPeptide1", action="store")
    MHChelixPeptide1 = file_name + "/visualisation/" + pdb_name + "_omitmap_MHChelixPeptide1.png"

    if ray:
        ray_tracer(MHChelixPeptide1, 1)
    else:
        ray_tracer(MHChelixPeptide1, 0)

    # Copy above scene but with a pink peptide
    pymol.cmd.color("magenta", "ps")
    pymol.cmd.util.cnc("ps")

    # Photo op here

    pymol.cmd.scene(key="MHChelixPeptide2", action="store")
    MHChelixPeptide2 = file_name + "/visualisation/" + pdb_name + "_omitmap_MHChelixPeptide2.png"

    if ray:
        ray_tracer(MHChelixPeptide2, 1)
    else:
        ray_tracer(MHChelixPeptide2, 0)
    # Save the session
    pymol.cmd.save(file_name + "/sessions/" + pdb_name + "_omitmap.pse")

    # Quit pymol
    pymol.cmd.quit()











