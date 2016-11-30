import argparse
import numpy as np
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
	Usage: python SoluteInCustomSolBox.py -x OCT.pdb -s UNK_9895E0.pdb -Nx OCT -Ns UNK -b 40 
	
	REQUIREMENTS:
	Preferably Anaconda python with following modules
	argparse
	numpy
	"""
    )
    parser.add_argument(
        "-s", "--pdb_solute", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument(
        "-x", "--pdb_solvent", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument(
        "-Ns", "--resid_solute", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument(
        "-Nx", "--resid_solvent", help="Submit PDB file from CHEMSPIDER or PubChem", type=str)
    parser.add_argument("-b", "--box_size", type=float,
                        help="SIZE of the CUBIC box in ANGSTROM")
    parser.add_argument("-nX", "--num_solv", type=int,
                        help="NUMBER of Molecules in CUBIC box")
    args = parser.parse_args()
    try:
        SOLVATE_SOLUTE(args.pdb_solvent,args.pdb_solute, args.resid_solvent,args.resid_solute, args.box_size,args.num_solv)
    except TypeError:
	print('For Help: python SoluteInCustomSolBox.py -h')
