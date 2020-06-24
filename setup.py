import setuptools
import os
import subprocess

#wget ccp4, untar it, move it into our library and activate the binaries
subprocess.run(["bash", "helpers/install_ccp4.sh"])

#add CCP4 files recursively when it comes to setup
def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

#todo see if there's a way to work out what libraries I can get away with deleting
ccp4_files = package_files('stacei_src/bin/ccp4-7.0')

#then install the package with setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="STACEI-WhalleyT",
    version="0.1",
    author="Tom Whalley",
    author_email="twhalley93@gmail.com",
    description="Python package for STACEI: structural analysis for TCR-pMHC complexes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/whalleyt/stacei",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU License",
        "Operating System :: Linux",
    ],
    entry_points={"console_scripts": [
        "STACEI = stacei_src.STACEI:main"
    ]},
    python_requires='>=3.6',
    install_requires=open("requirements.txt").read(),
    dependency_links=[
        "https://github.com/oxpig/ANARCI#egg=anarci",
        "https://github.com/schrodinger/pymol-open-source#egg=pymol"
    ],
    package_data={'stacei_src': ['bin/web/*', 'bin/data/*', 'bin/executables/*', 'bin/R/*'] + ccp4_files},
    zip_safe = False
)
