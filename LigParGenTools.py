from math import *
from random import randint, random
import numpy as np
import sys
import pandas as pd
############################################


def trans(coos, qi):
    rcoos = np.zeros((len(coos), len(qi)))
    for i in range(len(coos)):
        rcoos[i][0] = coos[i][0] + qi[0]
        rcoos[i][1] = coos[i][1] + qi[1]
        rcoos[i][2] = coos[i][2] + qi[2]
    return(rcoos)
#################### CENTER OF COOS #############################


def CenOfCoos(atoms, coos):
    com = np.zeros(3)
    for i in range(0, len(atoms)):
        com[0] += coos[i][0]
        com[1] += coos[i][1]
        com[2] += coos[i][2]
    ncom = [i / len(atoms) for i in com]
    rad = []
    for i in range(0, len(atoms)):
        coos[i][0] = coos[i][0] - ncom[0]
        coos[i][1] = coos[i][1] - ncom[1]
        coos[i][2] = coos[i][2] - ncom[2]
        rad.append(sqrt(coos[i][0]**2 + coos[i][1]**2 + coos[i][2]**2))
    return(coos, ncom, max(rad))
####################MOMENTS OF INERTIA #############################


def rotzyx(coos):
    phi = random() * 2.0 * np.pi
    csphi = np.cos(phi)
    snphi = np.sin(phi)
    chi = random() * 2.0 * np.pi
    cschi = np.cos(chi)
    snchi = np.sin(chi)
    cstta = (random() * 2) - 1
    sntta = np.sin(np.arccos(cstta))
    RXX = (cstta * csphi * cschi) - (snphi * snchi)
    RXY = -cstta * csphi * snchi - snphi * cschi
    RXZ = sntta * csphi
    RYX = cstta * snphi * cschi + csphi * snchi
    RYY = (-cstta * snphi * snchi) + (csphi * cschi)
    RYZ = sntta * snphi
    RZX = -sntta * cschi
    RZY = sntta * snchi
    RZZ = cstta
    rcoos = np.zeros((len(coos), 3))
    for i in range(len(coos)):
        rcoos[i][0] = coos[i][0] * RXX + coos[i][1] * RXY + coos[i][2] * RXZ
        rcoos[i][1] = coos[i][0] * RYX + coos[i][1] * RYY + coos[i][2] * RYZ
        rcoos[i][2] = coos[i][0] * RZX + coos[i][1] * RZY + coos[i][2] * RZZ
    return (rcoos)
########################## FOR READING BOXMAKER COOS ######################


def read_pdb(fname):
    fcon = open(fname).readlines()
    atoms = []
    coos = []
    pdb_lines = {}
    cons=[]
    for line in fcon:
        if ('ATOM' in line) or ('HETATM' in line):
            atoms.append(line[12:16].strip())
            coos.append(list(map(float, line[32:54].split())))
	if('CONECT' in line): cons.append(map(int,line.split()[1:]))
            # pdb_lines[line[12:16].strip()]=list(map(float,line[32:54].split()))
    return atoms, coos, cons

########################## FOR READING BOXMAKER COOS ######################


def pdb_lines(atoms, coos, ID, resid='UNK',offset=0):
    lines = []
    num = offset+1
    for (i, j) in zip(atoms, coos):
        lines.append('%-6s%5d %4s %3s A%4d    %8.3f%8.3f%8.3f\n' %
                     ('ATOM', num, i, resid, ID, j[0], j[1], j[2]))
        num += 1
#    lines.append('TER\n')
    return lines

def con_lines(cons,offset):
    lines = []
    num = offset
    for i in cons:
        lines.append('%-6s%5d%5d\n' %
                     ('CONECT',i[0]+num,i[1]+num))
    return lines

########################## FOR READING BOXMAKER COOS ######################


def BOX_MAKER(pdb_file, BOX_SIZE, res_name,NSolv=None):
    atMOL, csMOL,conMOL = read_pdb(pdb_file)
    csMOL, comMOL, radMOL = CenOfCoos(atMOL, csMOL)
#    BOX_SIZE = BOX_SIZE+ radMOL
    maxX = radMOL *2.0 
    maxY = radMOL *2.0
    maxZ = radMOL *2.0
    (nx, ny, nz) = (int(BOX_SIZE / maxX),
                    int(BOX_SIZE / maxY), int(BOX_SIZE / maxZ))
    gx = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=nx)
    gy = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=ny)
    gz = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=nz)

    BOX = {}
    ID = 1
    total_lines = []
    tot_con_lines = []
    for XC in gx:
        for YC in gy:
            for ZC in gz:
                csMOL = rotzyx(csMOL)
                csF = trans(csMOL, [XC, YC, ZC])
                BOX[ID] = {'ATS': atMOL, 'CS': csF}
                ID = ID + 1
    if NSolv: BOX = {I+1:BOX[I+1] for I in range(0,NSolv)}
    for ID in BOX.keys():
	total_lines += pdb_lines(BOX[ID]['ATS'], BOX[ID]['CS'], ID,resid=res_name[0:3],offset=(ID-1)*len(atMOL))
	tot_con_lines += con_lines(conMOL, offset=(ID-1)*len(atMOL))
    total_mols = len(BOX.keys())
    ofile = open('%d_%s_BOX.pdb'%(total_mols,res_name[0:3]), 'w+')
    ofile.write('LIGPARGEN GENERATED CUSTOM SOLVENT BOX \n')
    ofile.write("CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n"%(BOX_SIZE,BOX_SIZE,BOX_SIZE))
    for l in total_lines +['TER \n'] +tot_con_lines:
        ofile.write('%s' % l)
    ofile.write('END \n')
    ofile.close()
    return None
########################## FOR READING BOXMAKER COOS ######################
def SOLVATE_SOLUTE(X_file,S_file, X_name,S_name,BOX_SIZE,NSolv=None):
    '''
	X is SOLVENT and S is SOLUTE
    '''
    atMOL, csMOL,conMOL = read_pdb(X_file)
    csMOL, comMOL, radMOL = CenOfCoos(atMOL, csMOL)
    atS,     csS,  conS = read_pdb(S_file)
    csS,    comS,  radS = CenOfCoos(atS, csS)
    BOX_SIZE = BOX_SIZE+ radS
    maxX = radMOL*2.0 
    maxY = radMOL*2.0 
    maxZ = radMOL*2.0 
    (nx, ny, nz) = (int(BOX_SIZE / maxX),
                    int(BOX_SIZE / maxY), int(BOX_SIZE / maxZ))
    gx = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=nx)
    gy = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=ny)
    gz = np.linspace(-0.5 * BOX_SIZE, 0.5 * BOX_SIZE, num=nz)
    BOX = {}
    BOX[0]={'ATS': atS, 'CS': csS}
    tot_con_lines=[]
    resid_names={0:S_name[0:3]}
    resid_offset={0:0}
    ID = 1
    DIST_LIST={'ID':[0],'DIS':[0.00]}
    for XC in gx:
        for YC in gy:
            for ZC in gz:
		DIST_C = sqrt(XC**2+YC**2+ZC**2) 
		if DIST_C >= radS:
		    DIST_LIST['DIS'].append(DIST_C)
		    DIST_LIST['ID'].append(ID)
                    csMOL = rotzyx(csMOL)
                    csF = trans(csMOL, [XC, YC, ZC])
                    BOX[ID] = {'ATS': atMOL, 'CS': csF}
    		    resid_names[ID]=X_name[0:3]
    		    resid_offset[ID]=(ID-1)*len(atMOL)+len(atS)
                    ID = ID + 1
		else:
		    continue
    df_DL = pd.DataFrame(DIST_LIST)
    df_DL = df_DL.sort_values(['DIS'])
    total_mols = len(BOX.keys())
### ADDED THESE FOUR LINES TO ACCOUNT FOR NSOLV
    if NSolv and (total_mols > NSolv): 
	df_DL = df_DL.head(NSolv+1)
	BOX = {I:BOX[I] for I in df_DL.ID}
    else: 
	print "Increase the box size: \n BOX size is always bigger than final equilibrated box size"
	sys.exit()
    RID_list = BOX.keys() 
    RID_list.remove(0)
    total_lines=pdb_lines(BOX[0]['ATS'], BOX[0]['CS'], 1,resid=resid_names[0])
    tot_con_lines= con_lines(conS,0)
    for ID in RID_list:
	print ID
	total_lines += pdb_lines(BOX[ID]['ATS'], BOX[ID]['CS'], ID+1,resid=resid_names[ID],offset=resid_offset[ID])
	tot_con_lines += con_lines(conMOL, offset=resid_offset[ID])
    total_mols = len(BOX.keys())
    print total_mols
    ofile = open('%s_in_%s.pdb'%(S_name,X_name),'w+')
    ofile.write('LIGPARGEN GENERATED CUSTOM SOLVENT BOX \n')
    ofile.write("CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n"%(BOX_SIZE,BOX_SIZE,BOX_SIZE))
    for l in total_lines + ['TER \n']+ tot_con_lines:
        ofile.write('%s' % l)
    ofile.write('END \n')
    ofile.close()
    return None
########################## FOR READING BOXMAKER COOS ######################
