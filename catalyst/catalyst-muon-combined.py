# script-version: 2.0
# Combined Catalyst state: particle bunch + self-field mesh slices.
#
# OPALX registers one channel per particle container:
#   - "particle_bunch_0" (type "multimesh") — self-contained:
#       data/particles    — particle positions, IDs, attributes
#       data/boundingbox  — 2x2x2 bounding box
#       data/transform    — bunch-to-lab transform
#       data/mesh         — self-field PIC grid (E_self, rho, E_ext_*)
#
# This script visualises both: field slices through the mesh block +
# particle point gaussian coloured by B, with the accelerator geometry overlaid.
#
# Offline stand-ins:
#   - Particles+Fields: /home/biddisco/temp/muon/particle_bunch_0_snapshot/particle_bunch_0_snapshot.vtpc
#
# Under Catalyst the XMLMultiBlockDataReader (registrationName='particle_bunch_0')
# is replaced by the live Conduit channel automatically.
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
# Particle + mesh channel: particle_bunch_0
# =========================================================================
particle_bunch_0 = XMLMultiBlockDataReader(
    registrationName='particle_bunch_0',
    FileName=['/home/biddisco/temp/muon/particle_bunch_0_snapshot/particle_bunch_0_snapshot.vtpc'])

# Extract the mesh block (self-field PIC grid)
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

# =========================================================================
# Particle blocks: extract transform + particles
# =========================================================================

# Extract the transform block (Input 1 for the programmable filter)
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=particle_bunch_0)
extractBlock1.Set(
    Assembly='Hierarchy',
    Selectors=['//transform'],
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
extractBlock2 = ExtractBlock(registrationName='ExtractBlock2', Input=particle_bunch_0)
extractBlock2.Set(
    Assembly='Hierarchy',
    Selectors=['//particles'],
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
    ext_p = CreateExtractor('VTPC', particle_bunch_0,
                            registrationName='Snapshot_particle_bunch_0')
    ext_p.Trigger = 'Time Step'
    ext_p.Trigger.Frequency = _snapshot_step
    ext_p.Writer.FileName = 'particle_bunch_0_snapshot/{timestep:06d}.vtpc'

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
    e_self_LUT.RescaleTransferFunctionToDataRange()

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
