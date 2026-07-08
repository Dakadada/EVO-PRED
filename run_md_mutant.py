import sys
from openmm.app import *
from openmm import *
from openmm.unit import *
from mdtraj.reporters import XTCReporter

if len(sys.argv) < 4:
    print("Error: Missing arguments.")
    print("Usage: python run_md_mutant.py <prmtop> <inpcrd> <output_xtc>")
    sys.exit(1)

prmtop_file = sys.argv[1]
inpcrd_file = sys.argv[2]
output_xtc = sys.argv[3]

print(f"Loading topology: {prmtop_file}...")
prmtop = AmberPrmtopFile(prmtop_file)
print(f"Loading coordinates: {inpcrd_file}...")
inpcrd = AmberInpcrdFile(inpcrd_file)

# System setup matching standard production parameters
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
platform = Platform.getPlatformByName('CUDA')

simulation = Simulation(prmtop.topology, system, integrator, platform)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
simulation.context.setPositions(inpcrd.positions)

print("Executing energy minimization...")
simulation.minimizeEnergy()

# 5,000,000 steps at 2 fs time step = 10 ns trajectory
steps = 5000000
simulation.reporters.append(StateDataReporter(sys.stdout, 50000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(XTCReporter(output_xtc, 50000))

print(f"Starting {steps} steps (10 ns) of production MD...")
simulation.step(steps)
print("Production run complete.")