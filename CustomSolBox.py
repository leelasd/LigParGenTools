import argparse
import numpy as np
import os 
from LigParGenTools import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Converter.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
	SCRIPT TO CREATE CUSTOM SOLVENT BOXES FOR 
	OPENMM AND NAMD FROM LIGPARGEN FILES
	Created on Mon Nov 14 15:10:05 2016
	@author: Leela S. Dodda leela.dodda@yale.edu
	@author: William L. Jorgensen Lab 
	
	if using PDB file 
	Usage: python CustomSolBox.py -p OCT.pdb -b 45 -r OCT 
	Usage: python CustomSolBox.py -p CYH.pdb -r CYH -b 28 
	
	REQUIREMENTS:
	Preferably Anaconda python with following modules
	argparse
	numpy
	"""
    )
    parser.add_argument(
        "-p", "--pdb", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument(
        "-r", "--sol_name", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument("-b", "--box_size", type=float,
                        help="SIZE of the CUBIC box in ANGSTROM")
    parser.add_argument("-ns", "--num_solv", type=int,
                        help="NUMBER of Molecules in CUBIC box")
    args = parser.parse_args()
    try:
        BOX_MAKER(args.pdb, args.box_size,args.sol_name, args.num_solv)
    except TypeError:
	print('For Help: python CustomSolBox.py -h')
