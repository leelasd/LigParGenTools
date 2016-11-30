from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from dcdreporter import DCDReporter
import numpy as np

def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system

def Minimize(simulation,iters=0):
    simulation.minimizeEnergy(maxIterations=iters)
    position = simulation.context.getState(getPositions=True).getPositions()
    energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    PDBFile.writeFile(simulation.topology, position,
                          open('gasmin.pdb', 'w'))
    print 'Energy at Minima is %3.3f kcal/mol' % (energy._value * KcalPerKJ)
    return simulation

pdb = PDBFile('125_OCT_BOX.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('OCT.xml')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer)
system = OPLS_LJ(system)
## FOR NPT
TEMP = 300*kelvin
system.addForce(MonteCarloBarostat(1*bar,TEMP))
integrator = LangevinIntegrator(TEMP, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
#print('MINIMIZATION STARTED')
#simulation = Minimize(simulation,300)
#print('MINIMIZATION DONE')
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True,density=True))
#simulation.reporters.append(DCDReporter('myfile.dcd', 100, enforcePeriodicBox=False))
simulation.step(50000)
np_equ_pos = simulation.context.getState(getPositions=True).getPositions() 
PDBFile.writeFile(simulation.topology, np_equ_pos,open('NVT_EQ_FINAL.pdb', 'w'))
