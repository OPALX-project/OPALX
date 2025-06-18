## Simple Multithread test for OPALX
written by Jan Luca Binder

### Abstract
This test will generate the slurm script and the run OPALX with different amount of nodes and threads. As input it uses a give template input file and parameters that can be set in the script

### Generating the data
In [generate.py](generate.py) the amount of nodes and threads can be set. For this one has to use the `parameters variable`. The first argument decides the amount of nodes and the second argument of each array is the amount of threads
Example
```python
parameters = [
    # nodes, threads
    [ 1    , 1],
    [ 1    , 2],
    [ 1    , 4],
    [ 1    , 8],
    [ 1    , 16],
    [ 1    , 32],
]
```

Additionally in the begining of the script the following variables can be changed:
- `ROOT_FOLDER`: Decides where the script will generate the folders which include the template input file and slurm script
- `SLURM_FILE_NAME`: Name of the slurm script file that will be generated
- `INPUT_FILE`: File name of the input file. Note that the input file should be in the same directory as the script so it can copy it succesfully into the generated folders
- `OPALX_EXECUTABLE_FILE`: File path to the OPALX-executable


Running the script will then automatically generated all folders and slurm scripts. After which it will copy the input file to all different generated folders and sbatch the slurm files.

### Plotting
Running [test.py](test.py) will plot all times and errors of all generated data. Here the most important variables are
- `failed_col`: This will automatically be filled with all colounms where the error is bigger than `EPSILON`. Add your own coloumns and it will plot it aswell
- `time_col`: This decides which cols will be plotted in the times plot

Note that `test.py` should be run after the OPALX finished running

### Cleaning
Running `clean.sh` will remove all folders inside the out folder and remove all plots