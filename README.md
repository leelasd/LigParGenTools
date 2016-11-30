# LigParGenTools
Scripts to set up condensed phase simulations using of LigParGen Server 
----


Most of the FF such as AMBER, OPLS-AA and CHARMM are built using bottom up approach i.e, they are parameterized to reproduce condensed phase properties of organic building blocks of macromolecules. To validate the performance of these force fields one need to look at their performance in reproducing pure liquid properties. Creation of custom solvent boxes of desired size with desired number of solvent molecules is the first step in doing pure liquid simulations. This script is created for this purpose. 

## 1. Steps to setup Pure liquid simulations
 
 1. Choose NAMD, Gromacs or OpenMM for your simulation 
 2. Get parameter and topology files for solvent molecules (1-Octanol) using LigParGen Server
 3. Create a box of Octanol using `CustomSolBox.py` using the following command.
 
 ```
 python CustomSolBox.py -p OCT.pdb -b 65 -r OCT
 ```
 This creates a box of 125 Octanol molecules randomly oriented and needs to be minimized and equilibrated until density matches the experimental value. 

<img style="float: right;" src="https://github.com/leelasd/LigParGenTools/blob/master/Pliq_BOX.jpg" width="300" height="300" />

 4. Minimize the box and do  NPT Equibration for 2 nano seconds to get a good box. 
 
 ``` bash 
 python PLIQ_OPENMM.py 
 ```
