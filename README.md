# LigParGenTools
Scripts to set up condensed phase simulations using of LigParGen Server 
----

## Creating Custom Solvent Box For Pure liquid simulations
Most of the FF such as AMBER, OPLS-AA and CHARMM are built using bottom up approach i.e, they are parameterized to reproduce condensed phase properties of organic building blocks of macromolecules. To validate the performance of these force fields one need to look at their performance in reproducing pure liquid properties. Creation of custom solvent boxes of desired size with desired number of solvent molecules is the first step in doing pure liquid simulations. This script is created for this purpose. 
 
 ### Steps to setup Pure liquid simulation: 
 1. Choose NAMD, Gromacs or OpenMM for your simulation 
 2. Get parameter and topology files for solvent molecules (1-Octanol) using LigParGen Server
 3. Create a box of Octanol using `CustomSolBox.py` using the following command.
 
 ```
 python CustomSolBox.py -p OCT.pdb -b 65 -r OCT
 ```
 This creates a box of 125 Octanol molecules randomly oriented and needs to be minimized and equilibrated until density matches the experimental value. 
