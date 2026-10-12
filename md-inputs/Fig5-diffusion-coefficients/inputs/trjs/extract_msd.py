#!/usr/bin/env python3

import os.path as op
import pickle
import sys

import mdtraj as mdt
import numpy as np

sys.path.append(".")
from msd import (
    calculate_msd,
    calculate_msd_original,
    calculate_total_msd_original,
    calculate_total_msd_z_resolved,
    count_z_resolved,
)


def calculate_com_pos_bmi(positions, idcs, masses):
    assert idcs.shape[-1] == 3
    pos_bmi = positions[:, idcs.ravel(), :].reshape(
        positions.shape[0], -1, len(masses), 3
    )
    pos_com = (pos_bmi * masses[None, None, :, None]).sum(axis=2) / masses.sum()
    return pos_com


# topology and trj parameters
path_pdb = "pos.pdb"
assert op.isfile(path_pdb), path_pdb
top = mdt.load(path_pdb)
n_atoms = top.topology.n_atoms
idx_ele = top.topology.select("name C")
assert len(idx_ele) == 0
masses_bmi = np.array([67.07, 15.04, 57.12])
mass_pf6 = 144.960007
n_ion = n_atoms // 4
assert n_ion == 320
idx_pf6 = np.arange(n_ion * 3, n_atoms, 1)
i1 = np.arange(0, 3 * n_ion, 3)
i2 = np.arange(1, 3 * n_ion, 3)
i3 = np.arange(2, 3 * n_ion, 3)
idx_bmi = np.stack([i1, i2, i3]).T
# MSD parameters
stride = 4  # steps
zrange = (0.67, 11.82)  # z-bin range
zbins = 60  # z-bin number

# in and out directories from arguments
assert len(sys.argv) == 3
path_raw = sys.argv[1]
assert op.isdir(path_raw), path_raw
path_extracted = sys.argv[2]
assert op.isdir(path_extracted), path_extracted

""" read trajectory from lammps dump """
path_pickle = op.join(path_raw, "dump.pkl")
if op.isfile(path_pickle):
    print(f"Reading trj from pickle {path_pickle}")
    with open(path_pickle, "rb") as f:
        trj = pickle.load(f)
else:
    path_lammps = op.join(path_raw, "dump.lammpstrj")
    assert op.isfile(path_lammps), path_lammps
    print(f"Reading {path_lammps}")
    try:
        trj = mdt.load(path_lammps, top=top.topology)
    except:
        print("Failed to load trj wtih mdt")
        sys.exit(1)
    print(f"Writing trj to pickle {path_pickle}")
    with open(path_pickle, "wb") as f:
        pickle.dump(trj, f)

print("trj shape:", trj.xyz.shape)

""" Subtract COM """
idx_liquid = np.arange(trj.n_atoms - len(idx_ele))  # TODO
masses = np.array((list(masses_bmi) + [mass_pf6]) * len(idx_pf6))
total_mass = masses.sum()
com = (trj.xyz[:, idx_liquid, :] * masses[None, :, None]).sum(axis=1) / total_mass
# and subtract the first frame so that we get a difference of com position from 0-time
com -= (trj[0].xyz[:, idx_liquid, :] * masses[:, None]).sum(axis=1) / total_mass
# and correct the whole trajectories for the drift
trj.xyz[:, idx_liquid, :] -= com[:, None, :]

""" extract MSDs from trj """
path_msd = op.join(path_extracted, "msd.pickle")
if op.isfile(path_msd):
    print(f"MSDs already pickled in {path_msd}")
    # print(f"Reading MSDs from pickle {path_msd}")
    # with open(path_msd, "rb") as f:
    # msd_anion, msd_anion_original, msd_cation, msd_cation_original = pickle.load(f)
else:
    print("Anion MSD")
    # calculate z_resolved MSDs
    msd_anion = calculate_total_msd_z_resolved(
        trj.xyz, indices=idx_pf6, stride_origin=stride, zrange=zrange, zbins=zbins
    )
    # calculate MSDs using the original implementation by Ondrej for benchmarking
    msd_anion_original = calculate_total_msd_original(
        trj.xyz, indices=idx_pf6, stride_origin=stride
    )
    print("Cation MSD")

    # calculate z_resolved MSDs
    msd_cation = calculate_total_msd_z_resolved(
        calculate_com_pos_bmi(trj.xyz, idx_bmi, masses_bmi),
        indices=np.arange(len(idx_bmi)),
        stride_origin=stride,
        zrange=zrange,
        zbins=zbins,
    )
    # calculate MSDs using the original implementation by Ondrej for benchmarking
    msd_cation_original = calculate_total_msd_original(
        calculate_com_pos_bmi(trj.xyz, idx_bmi, masses_bmi),
        indices=np.arange(len(idx_bmi)),
        stride_origin=stride,
    )
    print(f"Writing MSDs to pickle {path_msd}")
    with open(path_msd, "wb") as f:
        pickle.dump((msd_anion, msd_anion_original, msd_cation, msd_cation_original), f)


path_count = op.join(path_extracted, "count.pickle")
if False and op.isfile(path_count):
    print(f"Counts already pickled in {path_msd}")
else:
    print(f"Writing counts to pickle {path_count}")
    counts_anion = count_z_resolved(trj.xyz, idx_pf6, zrange, zbins)
    counts_cation_atoms = count_z_resolved(trj.xyz, idx_bmi, zrange, zbins)
    counts_cation = count_z_resolved(
        calculate_com_pos_bmi(trj.xyz, idx_bmi, masses_bmi),
        np.arange(len(idx_bmi)),
        zrange,
        zbins,
    )
    # print(calculate_com_pos_bmi(trj.xyz, idx_bmi, masses_bmi))
    # print(counts_anion)
    # print(counts_cation)
    # print(counts_cation_atoms)
    with open(path_count, "wb") as f:
        pickle.dump((counts_anion, counts_cation, counts_cation_atoms), f)
