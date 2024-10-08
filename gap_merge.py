from simtk.openmm.app import PDBFile, ForceField, Simulation
from simtk.openmm import LangevinIntegrator, CustomExternalForce
from simtk.unit import kelvin, picoseconds, nanometers, kilojoule_per_mole, nanometer
import mdtraj as md
import numpy as np

gap_start = 23
gap_end = 30

# 1. Load the original structure and AlphaFold predicted structure
original_pdb = PDBFile('original_structure.pdb')
alphafold_pdb = PDBFile('alphafold_structure.pdb')

# 2. Structure alignment using MDTraj
original_traj = md.load_pdb('original_structure.pdb')
alphafold_traj = md.load_pdb('alphafold_structure.pdb')

# Superposition for alignment
alphafold_traj.superpose(original_traj, atom_indices=original_traj.topology.select('backbone'))

# Example: Using residues near the gap (from 20 to 35)
gap_nearby_residues = range(gap_start-10, gap_end+10)

# Select backbone atoms near the gap
gap_nearby_atoms = original_traj.topology.select(
    f'residue {gap_nearby_residues.start} to {gap_nearby_residues.stop} and backbone'
)

# Superposition using only atoms near the gap
alphafold_traj.superpose(original_traj, atom_indices=gap_nearby_atoms)

# 3. Replace the coordinates for missing residues
# Example: The missing residues are numbered from 23 to 30
missing_residues = range(gap_start, gap_end)

for residue in missing_residues:
    # Get the atoms for the corresponding residue in the original structure
    original_atoms = original_traj.topology.select(f'residue {residue}')
    # Replace the coordinates with those from the AlphaFold structure
    original_traj.xyz[0, original_atoms, :] = alphafold_traj.xyz[0, original_atoms, :]

# Get the new positions after replacement
new_positions = original_traj.openmm_positions(0)

# 4. Perform energy minimization and local MD simulation with OpenMM

# Load the force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Build the system
system = forcefield.createSystem(original_pdb.topology, constraints=HBonds)

# Set up constraints to only move atoms near the gap
# Compute the center of the gap
gap_center = np.mean([new_positions[i] for i in original_traj.topology.select(f'residue {residue}')], axis=0)

# Calculate distances from each atom to the gap center and constrain atoms further away
distance_threshold = 1.0 * nanometers  # 10Å
force_constant = 1000.0 * kilojoule_per_mole / nanometer**2
position_restraint = CustomExternalForce('0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)')
position_restraint.addPerParticleParameter('x0')
position_restraint.addPerParticleParameter('y0')
position_restraint.addPerParticleParameter('z0')
position_restraint.addGlobalParameter('k', force_constant)

for atom in original_pdb.topology.atoms():
    atom_position = new_positions[atom.index]
    distance = np.linalg.norm(atom_position - gap_center)
    if distance > distance_threshold:
        position_restraint.addParticle(atom.index, atom_position.value_in_unit(nanometers))

system.addForce(position_restraint)

# Set up the simulation
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(original_pdb.topology, system, integrator)
simulation.context.setPositions(new_positions)

# Perform energy minimization
simulation.minimizeEnergy()

# Run short MD simulation if needed
simulation.step(1000)

# Save the results
state = simulation.context.getState(getPositions=True)
with open('refined_structure.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
