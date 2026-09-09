#!/usr/bin/env python
"""
ParaView Programmable Filter for transforming either accelerator geometry or
particles between the bunch reference frame and the lab/world frame.

The transform data is expected to live on a single-point polydata block that
carries the point-data arrays:

  BunchOrigin, BunchRotRow0/1/2       : lab -> bunch transform
  BunchInvOrigin, BunchInvRotRow0/1/2 : bunch -> lab transform

That block can be provided in one of two ways:

  1. As a separate source (e.g. ippl_particles_electrons.transform).  Set it
     as Input 1.
  2. As the third block inside the particle multimesh collection itself.  In
     that case pass the whole collection as Input 1 and the filter will locate
     the block containing BunchOrigin automatically.

Input 0 is always the dataset you want to transform (geometry or particles).
"""

from paraview.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np

# ---------------------------------------------------------------------------
# User-tunable flag
# ---------------------------------------------------------------------------
# False: lab -> bunch (transform geometry, bunch stays centred)
# True:  bunch -> lab (transform particles to world coords)
#
# INVERSE is normally passed in by the caller: the loader in a Catalyst/ParaView
# state script (see catalyst/catalyst-muon-test-1.py) exec's this file with the
# parameter set. When this script is pasted directly into a Programmable Filter
# instead, set the flag here or export BUNCH_TRANSFORM_INVERSE=0|1.
import os
try:
    INVERSE
except NameError:
    INVERSE = os.environ.get('BUNCH_TRANSFORM_INVERSE', '1') not in ('0', 'false', 'False')

# ---------------------------------------------------------------------------
# Debug logging
# ---------------------------------------------------------------------------
_DBG_FILE = "/tmp/bunch_transform_filter_debug.txt"

def _dbg(msg):
    with open(_DBG_FILE, "a") as f:
        f.write(msg + "\n")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _type_name(obj):
    return type(obj).__name__


def _list_arrays(dataset, indent="  "):
    out = []
    if dataset is None or not hasattr(dataset, 'GetPointData'):
        return [indent + "(no GetPointData)"]
    pda = dataset.GetPointData()
    n = pda.GetNumberOfArrays()
    if n == 0:
        return [indent + "(no point arrays)"]
    for i in range(n):
        arr = pda.GetArray(i)
        name = pda.GetArrayName(i) if arr else "?"
        out.append(f"{indent}array[{i}] = {name}")
    return out


def _inspect_input(idx, dataset, indent=""):
    _dbg(f"{indent}Input {idx}: {_type_name(dataset)}")
    if dataset is None:
        return
    if hasattr(dataset, 'GetNumberOfPoints'):
        _dbg(f"{indent}  points: {dataset.GetNumberOfPoints()}")
    if hasattr(dataset, 'GetNumberOfCells'):
        _dbg(f"{indent}  cells: {dataset.GetNumberOfCells()}")
    if hasattr(dataset, 'GetNumberOfPartitionedDataSets'):
        n = dataset.GetNumberOfPartitionedDataSets()
        _dbg(f"{indent}  partitioned datasets: {n}")
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
            _inspect_input(f"{idx}[{i}]", block, indent + "  ")
        return
    if hasattr(dataset, 'GetNumberOfPartitions'):
        n = dataset.GetNumberOfPartitions()
        _dbg(f"{indent}  partitions: {n}")
        for i in range(n):
            block = dataset.GetPartition(i)
            _inspect_input(f"{idx}[{i}]", block, indent + "  ")
        return
    if hasattr(dataset, 'GetNumberOfBlocks'):
        n = dataset.GetNumberOfBlocks()
        _dbg(f"{indent}  blocks: {n}")
        for i in range(n):
            block = dataset.GetBlock(i)
            _inspect_input(f"{idx}[{i}]", block, indent + "  ")
        return
    for line in _list_arrays(dataset, indent + "  "):
        _dbg(line)


def _has_transform_arrays(dataset):
    """Return True if the leaf dataset has the BunchOrigin/BunchRotRow* arrays."""
    if dataset is None:
        return False
    if not hasattr(dataset, 'GetPointData'):
        return False
    pda = dataset.GetPointData()
    if pda is None:
        return False
    return (pda.GetArray("BunchOrigin") is not None
            and pda.GetArray("BunchRotRow0") is not None)


def _find_transform_block(dataset):
    """Recursively search a composite/multiblock dataset for the transform block."""
    if hasattr(dataset, 'GetNumberOfPartitionedDataSets'):
        n = dataset.GetNumberOfPartitionedDataSets()
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
            found = _find_transform_block(block)
            if found is not None:
                return found
        return None

    if hasattr(dataset, 'GetNumberOfPartitions'):
        n = dataset.GetNumberOfPartitions()
        for i in range(n):
            block = dataset.GetPartition(i)
            found = _find_transform_block(block)
            if found is not None:
                return found
        return None

    if hasattr(dataset, 'GetNumberOfBlocks'):
        n = dataset.GetNumberOfBlocks()
        for i in range(n):
            block = dataset.GetBlock(i)
            found = _find_transform_block(block)
            if found is not None:
                return found
        return None

    if _has_transform_arrays(dataset):
        return dataset

    return None


def _read_transform(dataset):
    """Extract origin and rotation rows from the transform block."""
    tf = _find_transform_block(dataset)
    if tf is None:
        raise RuntimeError(
            "Could not find a transform block with BunchOrigin/BunchRotRow* arrays. "
            "Make sure Input 1 is the ippl_particles_electrons.transform source "
            "or the full particle collection."
        )

    pda = tf.GetPointData()
    if INVERSE:
        origin = vtk_to_numpy(pda.GetArray('BunchInvOrigin'))[0]
        r0 = vtk_to_numpy(pda.GetArray('BunchInvRotRow0'))[0]
        r1 = vtk_to_numpy(pda.GetArray('BunchInvRotRow1'))[0]
        r2 = vtk_to_numpy(pda.GetArray('BunchInvRotRow2'))[0]
    else:
        origin = vtk_to_numpy(pda.GetArray('BunchOrigin'))[0]
        r0 = vtk_to_numpy(pda.GetArray('BunchRotRow0'))[0]
        r1 = vtk_to_numpy(pda.GetArray('BunchRotRow1'))[0]
        r2 = vtk_to_numpy(pda.GetArray('BunchRotRow2'))[0]

    return origin, np.vstack([r0, r1, r2])


def _transform_leaf(dataset, origin, R, path=""):
    """Transform points on a leaf polydata. Returns number of points processed."""
    pts = vtk_to_numpy(dataset.GetPoints().GetData())
    new_pts = (pts - origin) @ R.T
    dataset.GetPoints().SetData(numpy_to_vtk(new_pts))
    _dbg(f"  transformed leaf {path}: {pts.shape[0]} points")
    return pts.shape[0]


def _transform_dataset(dataset, origin, R, path=""):
    """Recursively transform points in a composite or leaf dataset."""
    if hasattr(dataset, 'GetNumberOfPartitionedDataSets'):
        n = dataset.GetNumberOfPartitionedDataSets()
        _dbg(f"{path} PartitionedDataSetCollection ({n} datasets)")
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
            _transform_dataset(block, origin, R, f"{path}[{i}]")
        return

    if hasattr(dataset, 'GetNumberOfPartitions'):
        n = dataset.GetNumberOfPartitions()
        _dbg(f"{path} PartitionedDataSet ({n} partitions)")
        for i in range(n):
            block = dataset.GetPartition(i)
            _transform_dataset(block, origin, R, f"{path}[{i}]")
        return

    if hasattr(dataset, 'GetNumberOfBlocks'):
        n = dataset.GetNumberOfBlocks()
        _dbg(f"{path} MultiBlock ({n} blocks)")
        for i in range(n):
            block = dataset.GetBlock(i)
            _transform_dataset(block, origin, R, f"{path}[{i}]")
        return

    if hasattr(dataset, 'GetPoints'):
        _transform_leaf(dataset, origin, R, path)
        return

    _dbg(f"{path} unrecognized dataset type: {_type_name(dataset)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
#_dbg("=" * 60)
input0 = self.GetInputDataObject(0, 0)
input1 = self.GetInputDataObject(0, 1)
#_inspect_input(0, input0)
#_inspect_input(1, input1)

origin, R = _read_transform(input1)
#_dbg(f"Transform origin: {origin}")
#_dbg(f"Transform rotation rows:\n{R}")

output = self.GetOutput()
output.DeepCopy(input0)
_transform_dataset(output, origin, R, "out")
#_dbg("Output generated")
