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
INVERSE = False  # False: lab -> bunch (transform geometry, bunch stays centred)
                 # True:  bunch -> lab (transform particles to world coords)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
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
    # vtkPartitionedDataSetCollection / vtkMultiBlockDataSet
    if hasattr(dataset, 'GetNumberOfPartitionedDataSets'):
        n = dataset.GetNumberOfPartitionedDataSets()
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
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


# ---------------------------------------------------------------------------
# Read the transform from Input 1
# ---------------------------------------------------------------------------
transform_input = self.GetInputDataObject(0, 1)
origin, R = _read_transform(transform_input)

# ---------------------------------------------------------------------------
# Apply the transform to the points of Input 0
# ---------------------------------------------------------------------------
input0 = self.GetInputDataObject(0, 0)
pts = vtk_to_numpy(input0.GetPoints().GetData())

# x' = R * (x - origin)
new_pts = (pts - origin) @ R.T

# ---------------------------------------------------------------------------
# Produce output
# ---------------------------------------------------------------------------
output = self.GetOutput()
output.ShallowCopy(input0)
output.GetPoints().SetData(numpy_to_vtk(new_pts))
