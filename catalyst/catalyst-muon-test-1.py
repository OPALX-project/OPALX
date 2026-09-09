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

# create a new 'XML PolyData Reader'
mue4_analytical_elementsvtp = XMLPolyDataReader(registrationName='mue4_analytical_elements.vtp', FileName=['/home/biddisco/temp/muon/mue4_analytical_elements.vtp'])
mue4_analytical_elementsvtp.Set(
    CellArrayStatus=['element_type'],
    TimeArray='None',
)

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=[mue4_analytical_elementsvtp, extractippl_particles_electronstransform])
programmableFilter1.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    PythonPath='',
    Script="""#!/usr/bin/env python
\"\"\"
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
\"\"\"

from paraview.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np

# ---------------------------------------------------------------------------
# User-tunable flag
# ---------------------------------------------------------------------------
INVERSE = True  # False: lab -> bunch (transform geometry, bunch stays centred)
                 # True:  bunch -> lab (transform particles to world coords)

# ---------------------------------------------------------------------------
# Debug logging
# ---------------------------------------------------------------------------
_DBG_FILE = "/tmp/bunch_transform_filter_debug.txt"

def _dbg(msg):
    pass


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _type_name(obj):
    return type(obj).__name__


def _list_arrays(dataset, indent="  "):
    out = []
    if dataset is None or not hasattr(dataset, \'GetPointData\'):
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
    if hasattr(dataset, \'GetNumberOfPoints\'):
        _dbg(f"{indent}  points: {dataset.GetNumberOfPoints()}")
    if hasattr(dataset, \'GetNumberOfCells\'):
        _dbg(f"{indent}  cells: {dataset.GetNumberOfCells()}")
    if hasattr(dataset, \'GetNumberOfPartitionedDataSets\'):
        n = dataset.GetNumberOfPartitionedDataSets()
        _dbg(f"{indent}  partitioned datasets: {n}")
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
            _inspect_input(f"{idx}[{i}]", block, indent + "  ")
        return
    if hasattr(dataset, \'GetNumberOfPartitions\'):
        n = dataset.GetNumberOfPartitions()
        _dbg(f"{indent}  partitions: {n}")
        for i in range(n):
            block = dataset.GetPartition(i)
            _inspect_input(f"{idx}[{i}]", block, indent + "  ")
        return
    if hasattr(dataset, \'GetNumberOfBlocks\'):
        n = dataset.GetNumberOfBlocks()
        _dbg(f"{indent}  blocks: {n}")
        for i in range(n):
            block = dataset.GetBlock(i)
            _inspect_input(f"{idx}[{i}]", block, indent + "  ")
        return
    for line in _list_arrays(dataset, indent + "  "):
        _dbg(line)


def _has_transform_arrays(dataset):
    \"\"\"Return True if the leaf dataset has the BunchOrigin/BunchRotRow* arrays.\"\"\"
    if dataset is None:
        return False
    if not hasattr(dataset, \'GetPointData\'):
        return False
    pda = dataset.GetPointData()
    if pda is None:
        return False
    return (pda.GetArray("BunchOrigin") is not None
            and pda.GetArray("BunchRotRow0") is not None)


def _find_transform_block(dataset):
    \"\"\"Recursively search a composite/multiblock dataset for the transform block.\"\"\"
    if hasattr(dataset, \'GetNumberOfPartitionedDataSets\'):
        n = dataset.GetNumberOfPartitionedDataSets()
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
            found = _find_transform_block(block)
            if found is not None:
                return found
        return None

    if hasattr(dataset, \'GetNumberOfPartitions\'):
        n = dataset.GetNumberOfPartitions()
        for i in range(n):
            block = dataset.GetPartition(i)
            found = _find_transform_block(block)
            if found is not None:
                return found
        return None

    if hasattr(dataset, \'GetNumberOfBlocks\'):
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
    \"\"\"Extract origin and rotation rows from the transform block.\"\"\"
    tf = _find_transform_block(dataset)
    if tf is None:
        raise RuntimeError(
            "Could not find a transform block with BunchOrigin/BunchRotRow* arrays. "
            "Make sure Input 1 is the ippl_particles_electrons.transform source "
            "or the full particle collection."
        )

    pda = tf.GetPointData()
    if INVERSE:
        origin = vtk_to_numpy(pda.GetArray(\'BunchInvOrigin\'))[0]
        r0 = vtk_to_numpy(pda.GetArray(\'BunchInvRotRow0\'))[0]
        r1 = vtk_to_numpy(pda.GetArray(\'BunchInvRotRow1\'))[0]
        r2 = vtk_to_numpy(pda.GetArray(\'BunchInvRotRow2\'))[0]
    else:
        origin = vtk_to_numpy(pda.GetArray(\'BunchOrigin\'))[0]
        r0 = vtk_to_numpy(pda.GetArray(\'BunchRotRow0\'))[0]
        r1 = vtk_to_numpy(pda.GetArray(\'BunchRotRow1\'))[0]
        r2 = vtk_to_numpy(pda.GetArray(\'BunchRotRow2\'))[0]

    return origin, np.vstack([r0, r1, r2])


def _transform_leaf(dataset, origin, R, path=""):
    \"\"\"Transform points on a leaf polydata. Returns number of points processed.\"\"\"
    pts = vtk_to_numpy(dataset.GetPoints().GetData())
    new_pts = (pts - origin) @ R.T
    dataset.GetPoints().SetData(numpy_to_vtk(new_pts))
    _dbg(f"  transformed leaf {path}: {pts.shape[0]} points")
    return pts.shape[0]


def _transform_dataset(dataset, origin, R, path=""):
    \"\"\"Recursively transform points in a composite or leaf dataset.\"\"\"
    if hasattr(dataset, \'GetNumberOfPartitionedDataSets\'):
        n = dataset.GetNumberOfPartitionedDataSets()
        _dbg(f"{path} PartitionedDataSetCollection ({n} datasets)")
        for i in range(n):
            block = dataset.GetPartitionedDataSet(i)
            _transform_dataset(block, origin, R, f"{path}[{i}]")
        return

    if hasattr(dataset, \'GetNumberOfPartitions\'):
        n = dataset.GetNumberOfPartitions()
        _dbg(f"{path} PartitionedDataSet ({n} partitions)")
        for i in range(n):
            block = dataset.GetPartition(i)
            _transform_dataset(block, origin, R, f"{path}[{i}]")
        return

    if hasattr(dataset, \'GetNumberOfBlocks\'):
        n = dataset.GetNumberOfBlocks()
        _dbg(f"{path} MultiBlock ({n} blocks)")
        for i in range(n):
            block = dataset.GetBlock(i)
            _transform_dataset(block, origin, R, f"{path}[{i}]")
        return

    if hasattr(dataset, \'GetPoints\'):
        _transform_leaf(dataset, origin, R, path)
        return

    _dbg(f"{path} unrecognized dataset type: {_type_name(dataset)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
_dbg("=" * 60)
input0 = self.GetInputDataObject(0, 0)
input1 = self.GetInputDataObject(0, 1)
_inspect_input(0, input0)
_inspect_input(1, input1)

origin, R = _read_transform(input1)
_dbg(f"Transform origin: {origin}")
_dbg(f"Transform rotation rows:\
{R}")

output = self.GetOutput()
output.DeepCopy(input0)
_transform_dataset(output, origin, R, "out")
_dbg("Output generated")
""",
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from programmableFilter1
programmableFilter1Display = Show(programmableFilter1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'element_type'
element_typeLUT = GetColorTransferFunction('element_type')
element_typeLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=5.0,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
programmableFilter1Display.Set(
    Representation='Wireframe',
    ColorArrayName=['CELLS', 'element_type'],
    LookupTable=element_typeLUT,
    Opacity=0.25,
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

# get color legend/bar for element_typeLUT in view renderView1
element_typeLUTColorBar = GetScalarBar(element_typeLUT, renderView1)
element_typeLUTColorBar.Set(
    Title='element_type',
    ComponentTitle='',
)

# set color bar visibility
element_typeLUTColorBar.Visibility = 1

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
programmableFilter1Display.SetScalarBarVisibility(renderView1, False)

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

# get opacity transfer function/opacity map for 'element_type'
element_typePWF = GetOpacityTransferFunction('element_type')
element_typePWF.Set(
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

    if os.environ.get('IPPL_CATALYST_PNG', '1') == '1':
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

