# script-version: 2.0
# Catalyst state for visualising the OPALX field channels exported via IPPL.
#
# In the new mesh-centric channel model, all fields on the bunch mesh are
# published as a single Conduit channel named "bunch_mesh" (type "mesh")
# with multiple cell-data arrays: E_self, B_self, rho, phi, E_external,
# B_external.  This script binds to that single channel and selects
# arrays by name for visualization.
#
# Run in-situ with IPPL_CATALYST_PNG=OFF so this script is the sole consumer of
# the field channels; otherwise IPPL's per-channel png_ext_*field.py extractors
# grab the same channels and Catalyst cannot bind them twice.
import os

import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# Render view
# ----------------------------------------------------------------
renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[1065, 665],
    AxesGrid='Grid Axes 3D Actor',
    OrientationAxesVisibility=0,
    UseColorPaletteForBackground=0,
    Background=[0.0, 0.0, 0.0],
    CameraPosition=[0.0, 0.0, 0.05],
    CameraFocalPoint=[0.0, 0.0, 0.0],
    CameraViewUp=[0.0, 1.0, 0.0],
    CameraParallelScale=0.02,
)

# ----------------------------------------------------------------
# Mesh channel stand-in (substituted by live Conduit channel in-situ).
# All field arrays (E_self, B_self, rho, E_external, B_external, phi)
# are cell data on this single uniform grid source.
#
# For offline development, snapshot VTI files can be generated with
# OPALX_CATALYST_DUMP_SNAPSHOT=N using the field pipeline.
# ----------------------------------------------------------------
bunch_mesh = XMLImageDataReader(
    registrationName='bunch_mesh',
    FileName=['/home/biddisco/temp/muon/bunch_mesh_snapshot/bunch_mesh_snapshot_0_0.vti'])
# Select all field arrays by default (cell data)
bunch_mesh.CellArrayStatus = ['E_self', 'B_self', 'rho', 'E_external', 'B_external', 'RankID']

# MergeBlocks to flatten per-rank partitions into a single dataset
bunch_merged = MergeBlocks(registrationName='bunch_mesh.MergedBlocks',
                           Input=bunch_mesh)
bunch_merged.MergePartitionsOnly = 1

# ----------------------------------------------------------------
# External B field: horizontal slice coloured by magnitude
# ----------------------------------------------------------------
b_external_cd2pd = CellDatatoPointData(registrationName='BExternalCD2PD',
                                       Input=bunch_merged)
b_external_cd2pd.CellDataArraytoprocess = ['B_external']

b_external_slice = Slice(registrationName='BExternalSlice', Input=b_external_cd2pd)
b_external_slice.SliceType = 'Plane'
b_external_slice.SliceType.Normal = [0.0, 1.0, 0.0]
b_external_slice.SliceType.Origin = [0.0, 0.0, 0.0]

b_external_sliceDisplay = Show(b_external_slice, renderView1, 'GeometryRepresentation')
b_external_sliceDisplay.Representation = 'Surface'
b_external_sliceDisplay.ColorArrayName = ['POINTS', 'B_external']
b_external_sliceDisplay.Opacity = 0.5

b_external_LUT = GetColorTransferFunction('B_external')
b_external_LUT.Set(AutomaticRescaleRangeMode='Grow and update every timestep')
b_external_LUTColorBar = GetScalarBar(b_external_LUT, renderView1)
b_external_LUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='B_external',
    ComponentTitle='Magnitude',
    AllowOverlappingLabels=1,
)
b_external_LUTColorBar.Visibility = 1

# ----------------------------------------------------------------
# External E field: vertical slice coloured by magnitude
# ----------------------------------------------------------------
e_external_cd2pd = CellDatatoPointData(registrationName='EExternalCD2PD',
                                       Input=bunch_merged)
e_external_cd2pd.CellDataArraytoprocess = ['E_external']

e_external_slice = Slice(registrationName='EExternalSlice', Input=e_external_cd2pd)
e_external_slice.SliceType = 'Plane'
e_external_slice.SliceType.Normal = [1.0, 0.0, 0.0]
e_external_slice.SliceType.Origin = [0.0, 0.0, 0.0]

e_external_sliceDisplay = Show(e_external_slice, renderView1, 'GeometryRepresentation')
e_external_sliceDisplay.Representation = 'Surface'
e_external_sliceDisplay.ColorArrayName = ['POINTS', 'E_external']
e_external_sliceDisplay.Opacity = 0.5

e_external_LUT = GetColorTransferFunction('E_external')
e_external_LUT.Set(AutomaticRescaleRangeMode='Grow and update every timestep')
e_external_LUTColorBar = GetScalarBar(e_external_LUT, renderView1)
e_external_LUTColorBar.Set(
    WindowLocation='Lower Right Corner',
    Title='E_external',
    ComponentTitle='Magnitude',
    AllowOverlappingLabels=1,
)
e_external_LUTColorBar.Visibility = 1

# ----------------------------------------------------------------
# Snapshot extractor for offline script development
# ----------------------------------------------------------------
# OPALX_CATALYST_DUMP_SNAPSHOT=N creates a data extractor that writes on iteration N
# (and every Nth iteration after that). Set to 0 or unset to disable.
_snapshot_step = 0
try:
    _snapshot_step = int(os.environ.get('OPALX_CATALYST_DUMP_SNAPSHOT', '0'))
except ValueError:
    _snapshot_step = 0

if _snapshot_step > 0:
    ext = CreateExtractor('VTPD', bunch_mesh, registrationName='Snapshot_bunch_mesh')
    ext.Trigger = 'Time Step'
    ext.Trigger.Frequency = _snapshot_step
    ext.Writer.FileName = 'bunch_mesh_snapshot/{timestep:06d}.vtpd'

# ----------------------------------------------------------------
# Per-timestep behaviour
# ----------------------------------------------------------------
def catalyst_execute(info):
    # Keep the transfer functions pinned to the current data range.
    b_external_LUT.RescaleTransferFunctionToDataRange()
    e_external_LUT.RescaleTransferFunctionToDataRange()

# ----------------------------------------------------------------
# PNG extractor
# ----------------------------------------------------------------
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'Time Step'
pNG1.Trigger.Frequency = 5
pNG1.Writer.Set(
    FileName='Fields_{timestep:06d}{camera}.png',
    ImageResolution=[1065, 665],
    Format='PNG',
)

# ----------------------------------------------------------------
# Catalyst options
# ----------------------------------------------------------------
from paraview import catalyst
options = catalyst.Options()
options.EnableCatalystLive = 1
# Execute the pipeline every timestep so that OPALX_CATALYST_DUMP_SNAPSHOT=N
# can be interpreted as the simulation iteration rather than the pipeline
# invocation count. Per-output frequencies (PNG, snapshot extractors) still
# throttle their own writes.
options.GlobalTrigger.Frequency = 1
options.CatalystLiveTrigger.Frequency = 5

# ----------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    SaveExtractsUsingCatalystOptions(options)
