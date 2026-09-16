# script-version: 2.0
# Catalyst state for visualising the OPALX field channels exported via IPPL.
#
# Channels consumed (only available when OPALX_CATALYST_FIELDS=ON):
#   ippl_vField_E_self      - space-charge electric field (bunch mesh)
#   ippl_vField_B_self      - space-charge magnetic field (bunch mesh)
#   ippl_vField_E_external  - beamline element electric field sampled on mesh
#   ippl_vField_B_external  - beamline element magnetic field sampled on mesh
#   ippl_sField_rho         - charge density on the bunch mesh
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
# Field channel stand-ins (substituted by live Conduit channels in-situ)
# The *_snapshot files can be generated with OPALX_CATALYST_DUMP_SNAPSHOT=N
# using the particle pipeline, which registers the same channels.
# ----------------------------------------------------------------
ippl_vField_E_self = XMLImageDataReader(
    registrationName='ippl_vField_E_self',
    FileName=['/home/biddisco/temp/muon/ippl_vField_E_self_snapshot/ippl_vField_E_self_snapshot_0_0.vti'])
ippl_vField_B_self = XMLImageDataReader(
    registrationName='ippl_vField_B_self',
    FileName=['/home/biddisco/temp/muon/ippl_vField_B_self_snapshot/ippl_vField_B_self_snapshot_0_0.vti'])
ippl_vField_E_external = XMLImageDataReader(
    registrationName='ippl_vField_E_external',
    FileName=['/home/biddisco/temp/muon/ippl_vField_E_external_snapshot/ippl_vField_E_external_snapshot_0_0.vti'])
ippl_vField_B_external = XMLImageDataReader(
    registrationName='ippl_vField_B_external',
    FileName=['/home/biddisco/temp/muon/ippl_vField_B_external_snapshot/ippl_vField_B_external_snapshot_0_0.vti'])
ippl_sField_rho = XMLImageDataReader(
    registrationName='ippl_sField_rho',
    FileName=['/home/biddisco/temp/muon/ippl_sField_rho_snapshot/ippl_sField_rho_snapshot_0_0.vti'])

# ----------------------------------------------------------------
# External B field: horizontal slice coloured by magnitude
# ----------------------------------------------------------------
b_external_cd2pd = CellDatatoPointData(registrationName='BExternalCD2PD',
                                       Input=ippl_vField_B_external)
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
                                       Input=ippl_vField_E_external)
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
# Snapshot extractors for offline script development
# ----------------------------------------------------------------
# OPALX_CATALYST_DUMP_SNAPSHOT=N creates data extractors that write on iteration N
# (and every Nth iteration after that). Set to 0 or unset to disable.
_snapshot_step = 0
try:
    _snapshot_step = int(os.environ.get('OPALX_CATALYST_DUMP_SNAPSHOT', '0'))
except ValueError:
    _snapshot_step = 0

if _snapshot_step > 0:
    _snapshot_sources = [
        ('ippl_vField_E_self', ippl_vField_E_self),
        ('ippl_vField_B_self', ippl_vField_B_self),
        ('ippl_vField_E_external', ippl_vField_E_external),
        ('ippl_vField_B_external', ippl_vField_B_external),
        ('ippl_sField_rho', ippl_sField_rho),
    ]
    for name, src in _snapshot_sources:
        ext = CreateExtractor('VTPD', src, registrationName=f'Snapshot_{name}')
        ext.Trigger = 'Time Step'
        ext.Trigger.Frequency = _snapshot_step
        ext.Writer.FileName = f'{name}_snapshot/{{timestep:06d}}.vtpd'

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
