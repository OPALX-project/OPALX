# script-version: 2.0
# Catalyst state generated using paraview version 6.1.2
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
    ViewSize=[1145, 665],
    CenterOfRotation=[-0.006948784008934089, -0.003960216380611253, -0.00011194456114840823],
    CameraPosition=[-1.010088270216936, 0.4501731382255094, -1.2336669887458869],
    CameraFocalPoint=[-0.2132561884024458, 0.04602854840922185, 0.13478171088121466],
    CameraViewUp=[0.08820561969932411, 0.9680948213971375, 0.2345467659924196],
    UseColorPaletteForBackground=0,
    Background=[0.0, 0.0, 0.0],
)

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PV Trivial Producer'
ippl_particles_electrons = PVTrivialProducer(registrationName='ippl_particles_electrons')

# Extract only the particle bunch (block_main), not the helper box.
electrons_bunch = ExtractBlock(
    registrationName='electrons_bunch',
    Input=ippl_particles_electrons,
    Selectors=['//block_main']
)

# Extract the transform block from the particle multimesh collection.
# The transform is embedded as block_transform, not published as a separate channel.
extractippl_particles_electronstransform = ExtractBlock(
    registrationName='electrons_transform',
    Input=ippl_particles_electrons,
    Selectors=['//block_transform']
)

# create a new 'XML Multi Block Data Reader' for the element geometry
# (built by tools/element_positions/element_positions.py from the
# data/<base>_ElementPositions.json that opalx writes during setup)
mue4_analytical_elementsvtm = XMLMultiBlockDataReader(registrationName='mue4_analytical_elements', FileName=['/home/biddisco/temp/muon/data/mue4_analytical_ElementPositions.vtm'])
mue4_analytical_elementsvtm.Set(
    CellArrayStatus=['etype'],
)

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=[mue4_analytical_elementsvtm, extractippl_particles_electronstransform])
programmableFilter1.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    PythonPath='',
    Script="""# The bunch-transform filter body lives in catalyst/bunch_transform_filter.py
# in the repo (the maintained copy of what used to be embedded here). It is
# exec'd with 'self' and the parameters below passed in.

# Transform direction:
#   True:  bunch -> lab (transform particles to world coords)
#   False: lab -> bunch (transform geometry, bunch stays centred)
INVERSE = True

import os
_filter_script = os.environ.get(
    'BUNCH_TRANSFORM_FILTER',
    '/home/biddisco/src/opalx/master/catalyst/bunch_transform_filter.py')
_filter_ns = {'self': self, 'INVERSE': INVERSE}
with open(_filter_script) as _fh:
    exec(compile(_fh.read(), _filter_script, 'exec'), _filter_ns)
""",
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from programmableFilter1
programmableFilter1Display = Show(programmableFilter1, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'etype'
etypeTF2D = GetTransferFunction2D('etype')

# get color transfer function/color map for 'etype'
# Categorical coloring by the per-block etype strings (field data on each
# element block of the vtm). Settings verbatim from the ParaView trace.
etypeLUT = GetColorTransferFunction('etype')
etypeLUT.Set(
    InterpretValuesAsCategories=1,
    AnnotationsInitialized=1,
    TransferFunction2D=etypeTF2D,
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=5.0,
    ),
    ScalarRangeInitialized=1.0,
    Annotations=['DRIFT', 'DRIFT', 'MONITOR', 'MONITOR', 'MULTIPOLE', 'MULTIPOLE', 'SBEND', 'SBEND', 'SOLENOID', 'SOLENOID'],
    IndexedColors=[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
    IndexedOpacities=[1.0, 1.0, 1.0, 1.0, 1.0],
    EnableOpacityMapping=1,
)

# trace defaults for the display properties.
programmableFilter1Display.Set(
    Representation='Surface',
    ColorArrayName=['FIELD', 'etype'],
    LookupTable=etypeLUT,
    Opacity=0.20,
)

# show data from electrons_bunch (particles only, no helper box)
ippl_particles_electronsDisplay = Show(electrons_bunch, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'B'
bLUT = GetColorTransferFunction('B')
bLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=0.00,
        range_max=0.90,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
ippl_particles_electronsDisplay.Set(
    Representation='Point Gaussian',
    ColorArrayName=['POINTS', 'B'],
    LookupTable=bLUT,
    Assembly='Assembly',
    GaussianRadius=0.0021,
#    IntensityScalars='B',
#    RadiusScalars='B',
#    TypeScalars='B',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
ippl_particles_electronsDisplay.ScaleTransferFunction.Points = [-0.0473624234221366, 0.0, 0.5, 0.0, 0.05172962957261628, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
ippl_particles_electronsDisplay.OpacityTransferFunction.Points = [-0.0473624234221366, 0.0, 0.5, 0.0, 0.05172962957261628, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for etypeLUT in view renderView1
etypeLUTColorBar = GetScalarBar(etypeLUT, renderView1)
etypeLUTColorBar.Set(
    AutoOrient=0,
    Orientation='Horizontal',
    WindowLocation='Lower Left Corner',
    Title='etype',
    ComponentTitle='',
    AllowOverlappingLabels=1,
)

# set color bar visibility
etypeLUTColorBar.Visibility = 1

# get color legend/bar for bLUT in view renderView1
bLUTColorBar = GetScalarBar(bLUT, renderView1)
bLUTColorBar.Set(
    WindowLocation='Upper Left Corner',
    Title='B',
    ComponentTitle='Magnitude',
    LabelFormat='{:.4f}',
    RangeLabelFormat='{:.4f}',
    AutomaticLabelFormat=0,
    DrawTickLabels=0,
)
# Shrink the scalar bar vertically to about 75% of its default height while
# keeping it in the top-left corner. Try the modern length property first;
# fall back to Position2 for older ParaView builds.
try:
    bLUTColorBar.ScalarBarLength = 0.25
except Exception:
    try:
        bLUTColorBar.Position2 = [0.13, 0.225]
    except Exception:
        pass

# set color bar visibility
bLUTColorBar.Visibility = 1

# show color legend
programmableFilter1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
ippl_particles_electronsDisplay.SetScalarBarVisibility(renderView1, True)

# Add a small fixed-width time stamp in the bottom-right corner.
timeLabel = Text(registrationName='TimeLabel')
timeLabel.Text = 't = 0000000.0 ps'
timeLabelDisplay = Show(timeLabel, renderView1, 'TextSourceRepresentation')
timeLabelDisplay.Set(
    FontSize=10,
    FontFamily='Courier',  # fixed-width so the label doesn't jitter
    Color=[1.0, 1.0, 1.0],
    WindowLocation='Lower Right Corner',
)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'B'
bPWF = GetOpacityTransferFunction('B')
bPWF.Set(
    Points=[0.001, 0.0, 0.5, 0.0, 0.2, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'etype'
etypePWF = GetOpacityTransferFunction('etype')
etypePWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup camera links
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()

# when using paraview gui
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = 'TimeStep'   # execute live every Catalyst timestep

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)

# ----------------------------------------------------------------
# In-situ Catalyst execution: save a PNG every timestep
# ----------------------------------------------------------------
import os

def catalyst_execute(info):
    """Render and save a screenshot for the current Catalyst timestep."""
    # Ensure the pipeline is up to date for this timestep
    renderView1.Update()

    # Auto-rescale the B color map to the current data range and keep fixed
    # decimal labels on the scalar bar.
    try:
        bLUT.RescaleTransferFunctionToDataRange()
        bLUTColorBar.LabelFormat = '{:.4f}'
        bLUTColorBar.RangeLabelFormat = '{:.4f}'
        bLUTColorBar.AutomaticLabelFormat = 0
        bLUTColorBar.DrawTickLabels = 0
    except Exception:
        pass

    # Build a filename from the Catalyst cycle/time
    cycle = int(info.cycle) if hasattr(info, 'cycle') else 0
    sim_time = float(info.time) if hasattr(info, 'time') else 0.0
    filename = f"/home/biddisco/temp/muon/catalyst_viz_{cycle:07d}.png"

    # Update the bottom-right time stamp (fixed-width so it doesn't jitter).
    try:
        timeLabel.Text = f"t = {sim_time:09.1f} ps"
    except Exception:
        pass

    if os.environ.get('IPPL_CATALYST_PNG', '1').upper() in ('1', 'ON', 'TRUE', 'YES'):
        SaveScreenshot(
            filename,
            renderView1,
            ImageResolution=[1024, 768],
            FontScaling='Scale fonts proportionally',
            OverrideColorPalette='',
            StereoMode='No change',
            TransparentBackground=0,
            CompressionLevel='5'
        )

