#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio.PDB import PDBParser, PDBIO, Select
import argparse
import math

class ResidueSelect(Select):
    def __init__(self, residues_to_remove):
        self.residues_to_remove = residues_to_remove

    def accept_residue(self, residue):
        return residue not in self.residues_to_remove

def parse_arguments():
    parser = argparse.ArgumentParser(description='PDBファイル内の指定されたギャップでCα距離が近い残基を検出し、削除します。')
    parser.add_argument('input_pdb', help='入力PDBファイルのパス')
    parser.add_argument('output_pdb', help='出力PDBファイルのパス')
    parser.add_argument('-d1', '--distance_gap1', type=float, default=8.0, help='ギャップ1残基の場合のCα間の距離閾値（Å）。デフォルトは8.0 Å')
    parser.add_argument('-d2', '--distance_gap2', type=float, default=10.0, help='ギャップ2残基の場合のCα間の距離閾値（Å）。デフォルトは10.0 Å')
    parser.add_argument('-d3', '--distance_gap3', type=float, default=12.0, help='ギャップ3残基の場合のCα間の距離閾値（Å）。デフォルトは12.0 Å')
    return parser.parse_args()

def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def main():
    args = parse_arguments()
    input_pdb = args.input_pdb
    output_pdb = args.output_pdb

    # ギャップごとの距離閾値を辞書で定義
    gap_thresholds = {
        1: args.distance_gap1,
        2: args.distance_gap2,
        3: args.distance_gap3
    }

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('structure', input_pdb)
    except Exception as e:
        print(f"エラー: PDBファイルの解析に失敗しました。詳細: {e}")
        sys.exit(1)

    residues_to_remove = set()

    for model in structure:
        for chain in model:
            # 通常の残基のみをリストに格納
            residues = [res for res in chain if res.id[0] == ' ']
            residues_sorted = sorted(residues, key=lambda r: r.id[1])

            for i, res in enumerate(residues_sorted):
                res_num = res.id[1]
                # チェックするギャップの範囲（1〜3）
                for gap in range(1, 4):  # gap=1,2,3
                    j = i + gap + 1  # gap=1 → j=i+2, gap=2 → j=i+3, gap=3 → j=i+4
                    if j >= len(residues_sorted):
                        break
                    next_res = residues_sorted[j]
                    next_res_num = next_res.id[1]
                    actual_gap = next_res_num - res_num - 1
                    if actual_gap == gap:
                        # Cα原子が存在するか確認
                        if 'CA' in res and 'CA' in next_res:
                            coord1 = res['CA'].get_coord()
                            coord2 = next_res['CA'].get_coord()
                            distance = calculate_distance(coord1, coord2)
                            threshold = gap_thresholds.get(gap, 8.0)  # デフォルト閾値は8.0 Å
                            if distance < threshold:
                                print(f"警告: チェイン {chain.id} の残基 {res_num} と {next_res_num} (ギャップ {gap}) のCα距離が {distance:.2f} Å で閾値 {threshold} Å 以下です。これらの残基を削除します。")
                                residues_to_remove.add(res)
                                residues_to_remove.add(next_res)

    if residues_to_remove:
        print(f"\n総削除残基数: {len(residues_to_remove)}")
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb, ResidueSelect(residues_to_remove))
        print(f"修正されたPDBファイルが {output_pdb} に保存されました。")
    else:
        print("問題のある残基は検出されませんでした。")

if __name__ == "__main__":
    main()
