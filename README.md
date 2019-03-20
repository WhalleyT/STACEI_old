# STACEI
STACEI: STructural Analysis of TCR-pEptide-MHC Interactions

STACEI is a tool, primarily written in Python and R designed for the analysis and exploration of T-cell receptor (TCR) to peptide:major histocompatibility complex (pMHC) crystal structures. For more information about the specific functionality please see the publication referenced in the authors section. The simplified summary of STACEI is shown below:

![](https://github.com/WhalleyT/STACEI/blob/master/repo_files/workflow.png?raw=true "STACEI workflow")


## Getting Started

### Prerequisites

STACEI leverages several pre-existing tools for its analyses. As such, the following must already be installed before using it:

* [CCP4 suite](http://www.ccp4.ac.uk/)

* [Python 2.7](https://www.python.org/download/releases/2.7/)

* [R](https://www.r-project.org/)

* [PyMol](https://www.schrodinger.com/suites/pymol)

* [ANARCI](http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php)

CCP4 can be installed by following the instruction guidlines. It must also be added to the path in order for Python to recognise it when making calls to it. This can be done by either running `./start` from within the CCP4 directory, or by running `source /applications/ccp4-x.y.x/bin/ccp4.setup-sh` and adding it to your BASH path.

Python and R generally come pre-installed on most Linux and Mac distributions, but if not they can be downloaded from the above links, or in Ubuntu with:
`sudo apt-get install python2.7`
`sudo apt-get install r-base`

PyMol can also be downloaded by following the above link. Conda users can install it using `conda install -c schrodinger pymol` and Ubuntu users can download it with `sudo apt-get install pymol`.

ANARCI must be downloaded by following the link manually.


### Installing
To install, simply clone this repository or download it. Cloning it can be done with:

`git clone github.com/whalleyt/STACEI`

### Docker
should you not want to manually install the tool but want a local version of the tool you can download it as a docker container.

First you must download an image of the tool like so:
`docker pull twhalley93/stacei`

Then create a container (an instance) of your container:
`docker run -dit stacei`

Now there is an active container, pass your pdb files into it:
`docker cp <pdb_files> <container_id>:/`

before finally entering the container:
`docker exec -it <container_id>`

Here one can run the tool much like you would on your normal bash shell, for example:
`python STACEI.py -F <structure> -R`

## Authors
Manuscript in  progress

## License

This project is licensed under the GNU GPL License - see the [LICENSE](LICENSE) file for details.
