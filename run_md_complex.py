import sys
from openmm.app import *
from openmm import *
from openmm.unit import *
from mdtraj.reporters import XTCReporter

print("Loading topology and coordinates...")
prmtop = AmberPrmtopFile('holo_complex.prmtop')
inpcrd = AmberInpcrdFile('holo_complex.inpcrd')

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')

simulation = Simulation(prmtop.topology, system, integrator, platform)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
simulation.context.setPositions(inpcrd.positions)

print("Minimizing energy...")
simulation.minimizeEnergy()

steps = 5000000
simulation.reporters.append(StateDataReporter(sys.stdout, 50000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(XTCReporter("holo_trajectory.xtc", 50000))

print(f"Running {steps} steps (10 ns)...")
simulation.step(steps)
print("Simulation complete.")