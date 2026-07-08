import sys
from openmm.app import *
from openmm import *
from openmm.unit import *
from mdtraj.reporters import XTCReporter
from pdbfixer import PDBFixer

if len(sys.argv) < 2:
    print("Error: Target PDB file required as argument.")
    sys.exit(1)

pdb_file = sys.argv[1]
output_prefix = pdb_file.split('.')[0]

print("Fixing PDB structure...")
fixer = PDBFixer(filename=pdb_file)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(True)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

print("Writing cleaned PDB...")
PDBFile.writeFile(fixer.topology, fixer.positions, open(f"{output_prefix}_fixed.pdb", 'w'))

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = Modeller(fixer.topology, fixer.positions)
modeller.addSolvent(forcefield, padding=1.0*nanometers)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

platform = Platform.getPlatformByName('CUDA')
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

print("Minimizing energy...")
simulation.minimizeEnergy()

steps = 500000
simulation.reporters.append(StateDataReporter(sys.stdout, 10000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(XTCReporter(f"{output_prefix}_trajectory.xtc", 10000))

print(f"Running {steps} steps...")
simulation.step(steps)
print("Simulation complete.")