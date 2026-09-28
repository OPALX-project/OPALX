# script-version: 2.0
# Catalyst state generated using paraview version 6.1.2
import os

import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[1065, 665],
    AxesGrid='Grid Axes 3D Actor',
    OrientationAxesVisibility=0,
    CenterOfRotation=[0.0030920839110432897, 0.003948736118858261, -0.0005461920502175784],
    CameraPosition=[-1.7426509604863896, 0.6465715578031851, -1.5098108192902613],
    CameraFocalPoint=[-0.3378740213177323, 0.16441397672705024, 0.3483630395754754],
    CameraViewUp=[0.09309152855820636, 0.9785930078078131, 0.18354752076874903],
    UseColorPaletteForBackground=0,
    Background=[0.0, 0.0, 0.0],
)

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML MultiBlock Data Reader'
# Stand-in for the live 'particle_bunch_0' particle channel. Under Catalyst
# this reader is replaced by the in-situ Conduit channel with the same
# registrationName.
particle_bunch_0 = XMLMultiBlockDataReader(
    registrationName='particle_bunch_0',
    FileName=['/home/biddisco/temp/muon/particle_bunch_0_snapshot/particle_bunch_0_snapshot.vtpc'])

# create a new 'Extract Block'
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=particle_bunch_0)
extractBlock1.Set(
    Assembly='Hierarchy',
    Selectors=['//transform'],
)

# create a new 'XML MultiBlock Data Reader'
mue4_analytical_ElementPositionsvtm = XMLMultiBlockDataReader(registrationName='mue4_analytical_ElementPositions.vtm', FileName=['/home/biddisco/temp/muon/data/mue4_analytical_ElementPositions.vtm'])
mue4_analytical_ElementPositionsvtm.Set(
    AssemblySelectors=['/Root'],
    CellArrayStatus=['etype'],
    PointArrayStatus=['TCoords', 'Normals', 'TubeNormals'],
    TimeArray='None',
)

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime(registrationName='AnnotateTime1')
annotateTime1.Format = 'Time: {time:07.1f}'

# create a new 'Extract Block'
extractBlock2 = ExtractBlock(registrationName='ExtractBlock2', Input=particle_bunch_0)
extractBlock2.Set(
    Assembly='Hierarchy',
    Selectors=['//particles'],
)

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=[mue4_analytical_ElementPositionsvtm, extractBlock1])
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

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from programmableFilter1
programmableFilter1Display = Show(programmableFilter1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'element_type'
element_typeLUT = GetColorTransferFunction('element_type')
element_typeLUT.Set(
    InterpretValuesAsCategories=1,
    AnnotationsInitialized=1,
    NanColor=[0.5, 0.5, 0.5],
    Annotations=['DRIFT', 'DRIFT', 'MONITOR', 'MONITOR', 'MULTIPOLE', 'MULTIPOLE', 'SBEND', 'SBEND', 'SOLENOID', 'SOLENOID'],
    IndexedColors=[0.6500037908554077, 0.6500037908554077, 0.6500037908554077, 0.9000076055526733, 0.10000763088464737, 0.7000076174736023, 0.10000763088464737, 0.5499961972236633, 0.2500038146972656, 0.9000076055526733, 0.10000763088464737, 0.10000763088464737, 1.0, 0.5499961972236633, 0.0],
    IndexedOpacities=[1.0, 1.0, 1.0, 1.0, 1.0],
)

# trace defaults for the display properties.
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

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
programmableFilter1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
programmableFilter1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data from extractBlock2
extractBlock2Display = Show(extractBlock2, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'B'
bLUT = GetColorTransferFunction('B')
bLUT.Set(
    AutomaticRescaleRangeMode='Grow and update every timestep',
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=0.5,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
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

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
extractBlock2Display.ScaleTransferFunction.Points = [-0.0010852295478314294, 0.0, 0.5, 0.0, 0.0010451635207817868, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
extractBlock2Display.OpacityTransferFunction.Points = [-0.0010852295478314294, 0.0, 0.5, 0.0, 0.0010451635207817868, 1.0, 0.5, 0.0]

# show data from annotateTime1
annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')

# trace defaults for the display properties.
annotateTime1Display.Set(
    WindowLocation='Lower Right Corner',
    FontFamily='Courier',
)

# setup the color legend parameters for each legend in this view

# get color legend/bar for bLUT in view renderView1
bLUTColorBar = GetScalarBar(bLUT, renderView1)
bLUTColorBar.Set(
    WindowLocation='Upper Left Corner',
    Title='B',
    ComponentTitle='Magnitude',
    AllowOverlappingLabels=1,
)

# set color bar visibility
bLUTColorBar.Visibility = 1

# get color legend/bar for element_typeLUT in view renderView1
element_typeLUTColorBar = GetScalarBar(element_typeLUT, renderView1)
element_typeLUTColorBar.Set(
    WindowLocation='Lower Left Corner',
    Title='Element',
    ComponentTitle='',
    AllowOverlappingLabels=1,
)

# set color bar visibility
element_typeLUTColorBar.Visibility = 1

# show color legend
programmableFilter1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
extractBlock2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'B'
bPWF = GetOpacityTransferFunction('B')
bPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'element_type'
element_typePWF = GetOpacityTransferFunction('element_type')

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
pNG1.Trigger = 'Time Step'

# init the 'Time Step' selected for 'Trigger'
pNG1.Trigger.Frequency = 5

# init the 'PNG' selected for 'Writer'
pNG1.Writer.Set(
    FileName='RenderView1_{timestep:06d}{camera}.png',
    ImageResolution=[1065, 665],
    Format='PNG',
)

# ----------------------------------------------------------------
# setup camera links
# ----------------------------------------------------------------

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
    ext = CreateExtractor('VTPC', particle_bunch_0,
                          registrationName='Snapshot_particle_bunch_0')
    ext.Trigger = 'Time Step'
    ext.Trigger.Frequency = _snapshot_step
    ext.Writer.FileName = 'particle_bunch_0_snapshot/{timestep:06d}.vtpc'

# ----------------------------------------------------------------
# Per-timestep behaviour
# ----------------------------------------------------------------
def catalyst_execute(info):
    pass

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.EnableCatalystLive = 1

# Execute every timestep so OPALX_CATALYST_DUMP_SNAPSHOT is interpreted as the
# simulation iteration. Individual extractor frequencies still throttle writes.
options.GlobalTrigger.Frequency = 5

# init the 'Time Step' selected for 'CatalystLiveTrigger'
options.CatalystLiveTrigger.Frequency = 5

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
