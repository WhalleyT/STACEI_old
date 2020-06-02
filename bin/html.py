from distutils.dir_util import copy_tree
from shutil import rmtree

import os

"""
This module is used to generate our output into a HTML file
"""

def create_copy(directory_path):
    """
    Copy our master copy to the output directory
    """

    if os.path.exists(directory_path):
        rmtree(directory_path)

    copy_tree("web/", directory_path)

def rewrite_file(file, keywords, name):
    data = []

    with open(file) as f:
        for line in f:
            data.append(line)


    #now replace our lines
    outdata = []

    for line in data:
        for keyword in keywords:
            if keyword in line:
                line = line.replace(keyword, keyword + name + "_")
        outdata.append(line)

    replacement_file = open(file, "w")

    for outline in outdata:
        replacement_file.write(outline)


def make_html(name):
    create_copy(name)

    #rewrite the general file
    print("writing file for general visualisations")
    rewrite_file(name + "/gen_vis.html", ["sessions/"], name)
    #rewrite the crossing angle html
    print("writing file for crossing angle")
    rewrite_file(name + "/cross_ang.html", ["crossingAngle/", "sessions/"], name)
     #todo check the onAxis file

