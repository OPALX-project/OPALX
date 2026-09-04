#!/usr/bin/env python
"""
ParaView Programmable Filter for transforming either accelerator geometry or
particles between the bunch reference frame and the lab/world frame.

Usage:
  1. Add a `Programmable Filter` to your ParaView pipeline.
  2. Set the filter to have **two inputs**:
       Input 0 : the dataset you want to transform (geometry or particles)
       Input 1 : the `BunchTransformData_<label>` source produced by OPALX
  3. Paste this script into the filter's `Script` field.
  4. Set `INVERSE = False` to transform lab geometry into the bunch frame
     (camera tracks the bunch).  Set `INVERSE = True` to transform bunch-frame
     particles into lab/world coordinates.

The transform data channel is a single-point polydata set that carries the
following point-data arrays:

  BunchOrigin, BunchRotRow0/1/2       : lab -> bunch transform
  BunchInvOrigin, BunchInvRotRow0/1/2 : bunch -> lab transform

These arrays are constant across the single point; index [0] is sufficient.
"""

from paraview.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np

# ---------------------------------------------------------------------------
# User-tunable flag
# ---------------------------------------------------------------------------
INVERSE = False  # False: lab -> bunch (transform geometry); True: bunch -> lab (transform particles)

# ---------------------------------------------------------------------------
# Read the transform from Input 1 (BunchTransformData_<label>)
# ---------------------------------------------------------------------------
transform_input = self.GetInputDataObject(0, 1)
pd = transform_input.GetPointData()

if INVERSE:
    origin = vtk_to_numpy(pd.GetArray('BunchInvOrigin'))[0]
    r0 = vtk_to_numpy(pd.GetArray('BunchInvRotRow0'))[0]
    r1 = vtk_to_numpy(pd.GetArray('BunchInvRotRow1'))[0]
    r2 = vtk_to_numpy(pd.GetArray('BunchInvRotRow2'))[0]
else:
    origin = vtk_to_numpy(pd.GetArray('BunchOrigin'))[0]
    r0 = vtk_to_numpy(pd.GetArray('BunchRotRow0'))[0]
    r1 = vtk_to_numpy(pd.GetArray('BunchRotRow1'))[0]
    r2 = vtk_to_numpy(pd.GetArray('BunchRotRow2'))[0]

R = np.vstack([r0, r1, r2])

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
