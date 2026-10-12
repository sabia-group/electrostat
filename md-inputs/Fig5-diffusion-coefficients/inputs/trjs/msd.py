import numpy as np
from tqdm.notebook import tqdm


def calculate_msd_original(data, stride_origin):
    """Calculate MSD of 1D data.
    Originally implemented by Ondrej for JPC Lett. ab initio water paper w/ T. Markland.
    """

    msd = np.zeros(len(data), float)
    counts = np.zeros(len(data), int)

    for i in range(0, len(data) - 1, stride_origin):
        size = len(data) - i
        msd[:size] += (data[i] - data[i:]) ** 2
        counts[:size] += 1

    return msd, counts


def calculate_total_msd_original(positions, indices, stride_origin):
    """Calculate MSD for a given species in a trajectory."""

    msd = np.zeros((positions.shape[0], 3), float)  # keep all cartesian axes
    counts = np.zeros((positions.shape[0], 3), float)

    for idx in tqdm(indices):
        for c in range(3):
            msd_iter = calculate_msd_original(positions[:, idx, c], stride_origin)
            msd[..., c] += msd_iter[0]
            counts[..., c] += msd_iter[1]

    return msd / counts, counts


def calculate_msd(data, stride_origin, resolve_by_property, prop, prop_range, nbins):
    """Calculate MSD of 1D data, optionally resolved by a property (e.g., z-coordinate)."""

    msd = np.zeros((len(data), nbins), float)
    counts = np.zeros((len(data), nbins), int)

    bins = np.linspace(prop_range[0], prop_range[1], nbins, endpoint=False)

    for i in range(0, len(data) - 1, stride_origin):
        size = len(data) - i
        target_bin = (
            np.digitize(prop[i], bins) - 1
        )  # which bin should this contribution fall into
        msd[:size, target_bin] += (
            data[i] - data[i:]
        ) ** 2  # situate the contribution into the 2D array
        counts[
            :size, target_bin
        ] += 1  # add counts into the right place in the 2D array

    return msd, counts


def calculate_total_msd_z_resolved(positions, indices, stride_origin, zrange, zbins):
    """Calculate MSD for a given species in a trajectory."""

    msd = np.zeros((positions.shape[0], zbins, 3), float)  # keep all cartesian axes
    counts = np.zeros((positions.shape[0], zbins, 3), float)  # keep all cartesian axes

    for idx in tqdm(indices):
        z0 = positions[:, idx, 2]
        for c in range(3):
            msd_iter = calculate_msd(
                positions[:, idx, c],
                stride_origin,
                resolve_by_property=True,
                prop=z0,
                prop_range=zrange,
                nbins=zbins,
            )
            msd[..., c] += msd_iter[0]
            counts[..., c] += msd_iter[1]

    return msd / counts, counts


def count_z_resolved(positions, indices, zrange, zbins):
    counts = np.zeros(zbins, int)
    bins = np.linspace(zrange[0], zrange[1], zbins, endpoint=False)
    for idx in tqdm(indices):
        z0 = positions[:, idx, 2]
        for i in range(len(z0)):
            target_bin = np.digitize(z0[i], bins) - 1
            counts[target_bin] += 1
    return counts
