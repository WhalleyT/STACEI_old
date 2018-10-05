# STACEI
STACEI: STructural Analysis of TCR-pEptide-MHC Interactions

STACEI is a tool, primarily written in Python and R designed for the analysis and exploration of T-cell receptor (TCR) to peptide:major histocompatibility complex (pMHC) crystal structures. For more information about the specific functionality please see <citation>.


## Getting Started

### Prerequisites

STACEI leverages several pre-existing tools for its analyses. As such, the following must already be installed before using it:

* [CCP4 suite](http://www.ccp4.ac.uk/)

* [Python 2.7](https://www.python.org/download/releases/2.7/)

* [R](https://www.r-project.org/)

* [PyMol](https://www.schrodinger.com/suites/pymol)

* [ANARCI](http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php)

CCP4 can be installed by following the instruction guidlines. It must also be added to the path in order for Python to recognise it when making calls to it. This can be done by either running `./start` from within the CCP4 directory, or by running `source /applications/ccp4-x.y.x/bin/ccp4.setup-sh` and adding it to your BASH path.


### Installing
To install, simply clone this repository or download it. Cloning it can be done with:

`git clone github.com/whalleyt/STACEI` 

## Authors
Manuscript in  progress

## License

This project is licensed under the GNU GPL License - see the [LICENSE](LICENSE) file for details.

