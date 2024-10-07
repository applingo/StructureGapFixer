from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
import numpy as np


gap_start = 23
gap_end = 30

# 1. 元の構造とAlphaFoldで予測した構造を読み込む
original_pdb = PDBFile('original_structure.pdb')
alphafold_pdb = PDBFile('alphafold_structure.pdb')

# 2. 構造のアライメント（MDTrajを使用）
original_traj = md.load_pdb('original_structure.pdb')
alphafold_traj = md.load_pdb('alphafold_structure.pdb')

# スーパーポジションによるアライメント
alphafold_traj.superpose(original_traj, atom_indices=original_traj.topology.select('backbone'))

# 例として、ギャップの近傍残基として20から35を使用
gap_nearby_residues = range(gap_start-10, gap_end+10)

# ギャップ近傍のバックボーン原子を選択
gap_nearby_atoms = original_traj.topology.select(
    f'residue {gap_nearby_residues.start} to {gap_nearby_residues.stop} and backbone'
)

# スーパーポジションによるアライメント（ギャップ近傍の原子のみ）
alphafold_traj.superpose(original_traj, atom_indices=gap_nearby_atoms)



# 3. 欠損部分の残基番号を指定して座標を置換
# 例として、欠損部分が残基番号23から30とします
missing_residues = range(gap_start, gap_end)

for residue in missing_residues:
    # 元の構造の該当残基の原子を取得
    original_atoms = original_traj.topology.select(f'residue {residue}')
    # AlphaFold構造の対応する原子の座標で置換
    original_traj.xyz[0, original_atoms, :] = alphafold_traj.xyz[0, original_atoms, :]

# 置換後の座標を取得
new_positions = original_traj.openmm_positions(0)

# 4. OpenMMでエネルギー最小化と局所的なMDシミュレーションを行う

# フォースフィールドの読み込み
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# 系の構築
system = forcefield.createSystem(original_pdb.topology, constraints=HBonds)

# ギャップ周辺の原子のみを可動にするための拘束を設定
# ギャップ中心の座標を計算
gap_center = np.mean([new_positions[i] for i in original_traj.topology.select(f'residue {residue}')], axis=0)

# 各原子との距離を計算し、一定距離以上の原子に拘束をかける
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

# シミュレーションの設定
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(original_pdb.topology, system, integrator)
simulation.context.setPositions(new_positions)

# エネルギー最小化
simulation.minimizeEnergy()

# 必要に応じて短いMDシミュレーションを実行
simulation.step(1000)

# 結果の保存
state = simulation.context.getState(getPositions=True)
with open('refined_structure.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
