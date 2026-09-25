# script-version: 2.0
# Combined Catalyst state: particle bunch (MUON) + field slices (bunch_mesh).
#
# When OPALX_CATALYST_FIELDS=ON, OPALX registers two channels:
#   - "bunch_mesh"  (type "mesh")      — uniform grid with E_self, B_self, rho,
#                                        phi, E_external, B_external as cell data
#   - "MUON"        (type "multimesh") — particle block_main + block_transform
#
# This script visualises both: field slices through the mesh + particle point
# gaussian coloured by B, with the accelerator geometry overlaid.
#
# Offline stand-ins:
#   - Fields:  /home/biddisco/temp/muon/bunch_mesh_snapshot/bunch_mesh_snapshot_0_0.vti
#   - Particles: /home/biddisco/temp/muon/MUON_snapshot/MUON_snapshot.vtm
#
# Under Catalyst the XMLImageDataReader (registrationName='bunch_mesh') and
# XMLMultiBlockDataReader (registrationName='MUON') are replaced by the live
# Conduit channels automatically.
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
    CenterOfRotation=[0.0030920839110432897, 0.003948736118858581, -0.0005461920502175784],
    CameraPosition=[-1.7426509604863896, 0.6465715578031851, -1.5098108192902613],
    CameraFocalPoint=[-0.3378740213177323, 0.16441397672705024, 0.3483630395754754],
    CameraViewUp=[0.09309152855820636, 0.9785930078078131, 0.18354752076874903],
    UseColorPaletteForBackground=0,
    Background=[0.0, 0.0, 0.0],
)

# =========================================================================
# Field channel: bunch_mesh
# =========================================================================
bunch_mesh = XMLImageDataReader(
    registrationName='bunch_mesh',
    FileName=['/home/biddisco/temp/muon/bunch_mesh_snapshot/bunch_mesh_snapshot_0_0.vti'])
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

# =========================================================================
# Particle channel: MUON
# =========================================================================
muon_particles = XMLMultiBlockDataReader(
    registrationName='MUON',
    FileName=['/home/biddisco/temp/muon/MUON_snapshot/MUON_snapshot.vtm'])

# Extract the transform block (Input 1 for the programmable filter)
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=muon_particles)
extractBlock1.Set(
    Assembly='Hierarchy',
    Selectors=['//block_transform'],
)

# Accelerator geometry (unchanged file source, not a Catalyst channel)
mue4_analytical_ElementPositionsvtm = XMLMultiBlockDataReader(
    registrationName='mue4_analytical_ElementPositions.vtm',
    FileName=['/home/biddisco/temp/muon/data/mue4_analytical_ElementPositions.vtm'])
mue4_analytical_ElementPositionsvtm.Set(
    AssemblySelectors=['/Root'],
    CellArrayStatus=['etype'],
    PointArrayStatus=['TCoords', 'Normals', 'TubeNormals'],
    TimeArray='None',
)

# Annotate Time
annotateTime1 = AnnotateTime(registrationName='AnnotateTime1')
annotateTime1.Format = 'Time: {time:07.1f}'

# Extract the main particle block
extractBlock2 = ExtractBlock(registrationName='ExtractBlock2', Input=muon_particles)
extractBlock2.Set(
    Assembly='Hierarchy',
    Selectors=['//block_main'],
)

# Programmable Filter: transform particles to lab frame
programmableFilter1 = ProgrammableFilter(
    registrationName='ProgrammableFilter1',
    Input=[mue4_analytical_ElementPositionsvtm, extractBlock1])
programmableFilter1.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    CopyArrays=1,
    PythonPath='',
    Script="""INVERSE = True # parameter passed into the filter body
import os

_filter_script = os.environ.get(
'BUNCH_TRANSFORM_FILTER','/home/biddisco/src/opalx/insitu/catalyst/bunch_transform_filter.py')
_filter_ns = {'self': self, 'INVERSE': INVERSE}
with open(_filter_script) as _fh:
    exec(compile(_fh.read(), _filter_script, 'exec'), _filter_ns)""",
)

# Show transformed geometry (low opacity)
programmableFilter1Display = Show(programmableFilter1, renderView1, 'GeometryRepresentation')

element_typeLUT = GetColorTransferFunction('element_type')
element_typeLUT.Set(
    InterpretValuesAsCategories=1,
    AnnotationsInitialized=1,
    NanColor=[0.5, 0.5, 0.5],
    Annotations=['DRIFT', 'DRIFT', 'MONITOR', 'MONITOR', 'MULTIPOLE', 'MULTIPOLE', 'SBEND', 'SBEND', 'SOLENOID', 'SOLENOID'],
    IndexedColors=[0.6500037908554077, 0.6500037908554077, 0.6500037908554077, 0.9000076055526733, 0.10000763088464737, 0.7000076174736023, 0.10000763088464737, 0.5499961972236633, 0.2500038146972656, 0.9000076055526733, 0.10000763088464737, 0.10000763088464737, 1.0, 0.5499961972236633, 0.0],
    IndexedOpacities=[1.0, 1.0, 1.0, 1.0, 1.0],
)

programmableFilter1Display.Set(
    Representation='Surface',
    ColorArrayName=['FIELD', 'element_type'],
    LookupTable=element_typeLUT,
    Opacity=0.1,
    SelectNormalArray='TubeNormals',
    ComputePointNormals=1,
    Splitting=0,
    EdgeColor=[1.0, 1.0, 1.0],
    EdgeOpacity=0.0,
    Assembly='Hierarchy',
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
    IntensityScalars='Normals',
    RadiusScalars='Normals',
    TypeScalars='Normals',
)

# Show particles coloured by |B|
extractBlock2Display = Show(extractBlock2, renderView1, 'UnstructuredGridRepresentation')

bLUT = GetColorTransferFunction('B')
bLUT.Set(
    AutomaticRescaleRangeMode='Grow and update every timestep',
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=0.5,
    ),
    ScalarRangeInitialized=1.0,
)

extractBlock2Display.Set(
    Representation='Point Gaussian',
    ColorArrayName=['POINTS', 'B'],
    LookupTable=bLUT,
    Assembly='Hierarchy',
    GaussianRadius=0.0025,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
    IntensityScalars='B',
    RadiusScalars='B',
    TypeScalars='B',
)

bLUTColorBar = GetScalarBar(bLUT, renderView1)
bLUTColorBar.Set(
    WindowLocation='Upper Left Corner',
    Title='B',
    ComponentTitle='Magnitude',
    AllowOverlappingLabels=1,
)
bLUTColorBar.Visibility = 1

element_typeLUTColorBar = GetScalarBar(element_typeLUT, renderView1)
element_typeLUTColorBar.Set(
    WindowLocation='Lower Left Corner',
    Title='Element',
    ComponentTitle='',
    AllowOverlappingLabels=1,
)
element_typeLUTColorBar.Visibility = 1

programmableFilter1Display.SetScalarBarVisibility(renderView1, True)
extractBlock2Display.SetScalarBarVisibility(renderView1, True)

# Show time annotation
annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')
annotateTime1Display.Set(
    WindowLocation='Lower Right Corner',
    FontFamily='Courier',
)

# ----------------------------------------------------------------
# Opacity transfer functions
# ----------------------------------------------------------------
bPWF = GetOpacityTransferFunction('B')
bPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

element_typePWF = GetOpacityTransferFunction('element_type')

# =========================================================================
# Snapshot extractors for offline script development
# =========================================================================
_snapshot_step = 0
try:
    _snapshot_step = int(os.environ.get('OPALX_CATALYST_DUMP_SNAPSHOT', '0'))
except ValueError:
    _snapshot_step = 0

if _snapshot_step > 0:
    ext_p = CreateExtractor('VTPC', muon_particles,
                            registrationName='Snapshot_MUON')
    ext_p.Trigger = 'Time Step'
    ext_p.Trigger.Frequency = _snapshot_step
    ext_p.Writer.FileName = 'MUON_snapshot/{timestep:06d}.vtpc'

    ext_f = CreateExtractor('VTPD', bunch_mesh,
                            registrationName='Snapshot_bunch_mesh')
    ext_f.Trigger = 'Time Step'
    ext_f.Trigger.Frequency = _snapshot_step
    ext_f.Writer.FileName = 'bunch_mesh_snapshot/{timestep:06d}.vtpd'

# ----------------------------------------------------------------
# PNG extractor
# ----------------------------------------------------------------
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
pNG1.Trigger = 'Time Step'
pNG1.Trigger.Frequency = 5
pNG1.Writer.Set(
    FileName='Combined_{timestep:06d}{camera}.png',
    ImageResolution=[1065, 665],
    Format='PNG',
)

# ----------------------------------------------------------------
# Per-timestep behaviour
# ----------------------------------------------------------------
def catalyst_execute(info):
    # Rescale LUTs to current data range (fields move significantly)
    bLUT.RescaleTransferFunctionToDataRange()
    e_external_LUT.RescaleTransferFunctionToDataRange()

# ----------------------------------------------------------------
# Catalyst options
# ----------------------------------------------------------------
from paraview import catalyst
options = catalyst.Options()
options.EnableCatalystLive = 1

# Execute every timestep so snapshot extractors interpret N correctly.
options.GlobalTrigger.Frequency = 1
options.CatalystLiveTrigger.Frequency = 5

# ----------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    SaveExtractsUsingCatalystOptions(options)
