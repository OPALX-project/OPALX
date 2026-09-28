# script-version: 2.0
# Catalyst state for visualising the OPALX self-field mesh.
#
# The self-field PIC grid (E_self, rho, E_ext_*) is published as a "mesh"
# block inside the "particle_bunch_0" multimesh channel.  This script
# extracts the mesh block and renders field slices.
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
# Particle + mesh channel: particle_bunch_0
# The mesh block contains self-field arrays (E_self, rho, E_ext_*).
# ----------------------------------------------------------------
particle_bunch_0 = XMLMultiBlockDataReader(
    registrationName='particle_bunch_0',
    FileName=['/home/biddisco/temp/muon/particle_bunch_0_snapshot/particle_bunch_0_snapshot.vtpc'])

# Extract the mesh block
mesh_block = ExtractBlock(registrationName='MeshBlock', Input=particle_bunch_0)
mesh_block.Set(
    Assembly='Hierarchy',
    Selectors=['//mesh'],
)

# MergeBlocks to flatten per-rank partitions into a single dataset
mesh_merged = MergeBlocks(registrationName='mesh.MergedBlocks',
                          Input=mesh_block)
mesh_merged.MergePartitionsOnly = 1

# ----------------------------------------------------------------
# Self-field E: horizontal slice coloured by magnitude
# ----------------------------------------------------------------
e_self_cd2pd = CellDatatoPointData(registrationName='ESelfCD2PD',
                                   Input=mesh_merged)
e_self_cd2pd.CellDataArraytoprocess = ['E_self']

e_self_slice = Slice(registrationName='ESelfSlice', Input=e_self_cd2pd)
e_self_slice.SliceType = 'Plane'
e_self_slice.SliceType.Normal = [0.0, 1.0, 0.0]
e_self_slice.SliceType.Origin = [0.0, 0.0, 0.0]

e_self_sliceDisplay = Show(e_self_slice, renderView1, 'GeometryRepresentation')
e_self_sliceDisplay.Representation = 'Surface'
e_self_sliceDisplay.ColorArrayName = ['POINTS', 'E_self']
e_self_sliceDisplay.Opacity = 0.5

e_self_LUT = GetColorTransferFunction('E_self')
e_self_LUT.Set(AutomaticRescaleRangeMode='Grow and update every timestep')
e_self_LUTColorBar = GetScalarBar(e_self_LUT, renderView1)
e_self_LUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='E_self',
    ComponentTitle='Magnitude',
    AllowOverlappingLabels=1,
)
e_self_LUTColorBar.Visibility = 1

# ----------------------------------------------------------------
# Snapshot extractor for offline script development
# ----------------------------------------------------------------
_snapshot_step = 0
try:
    _snapshot_step = int(os.environ.get('OPALX_CATALYST_DUMP_SNAPSHOT', '0'))
except ValueError:
    _snapshot_step = 0

if _snapshot_step > 0:
    ext = CreateExtractor('VTPC', particle_bunch_0, registrationName='Snapshot_particle_bunch_0')
    ext.Trigger = 'Time Step'
    ext.Trigger.Frequency = _snapshot_step
    ext.Writer.FileName = 'particle_bunch_0_snapshot/{timestep:06d}.vtpc'

# ----------------------------------------------------------------
# Per-timestep behaviour
# ----------------------------------------------------------------
def catalyst_execute(info):
    e_self_LUT.RescaleTransferFunctionToDataRange()

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
options.GlobalTrigger.Frequency = 1
options.CatalystLiveTrigger.Frequency = 5

# ----------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    SaveExtractsUsingCatalystOptions(options)
