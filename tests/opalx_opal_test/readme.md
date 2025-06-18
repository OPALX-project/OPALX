## Opal - OpalX tests
written by Jan Luca Binder

### Abstract
This script will automatically create the input files and slurm script for OPALX and then run it for different particle amounts. Afterwhich it will compare the solution of the simulation with a given reference solution. 

### Introduction

To add a new test open [test.py](test.py). At the very top of the script you'll find 3 important global variables
- OPALX_EXECTUABLE_FILE
- parameters
- plotting_cols

### Executable path
First change the OPALX_EXECTUABLE_FILE string to the location of your the executable.
Example
```python
OPALX_EXECUTABLE_FILE = "/data/user/<user name>/opalx-elements/build/src/opalx"
```

### Parameters
Next you can edit the parameters. The name of the run is the key of each entry in the map. Each entry should have two further entries
1. An array called `amounts`: This decides how many particles will be used for each simulation
2. A file path called `ref`: The script excpects the reference stat files to be in a subfolder called `references`. Add your file there and add the name of the file to the map

Example
```python
parameters = {
    "run1" : {
        "amount" : ["1e2", "1e3", "1e4", "1e5"],
        "ref" : "ref-10steps.stat"
    },
    "run2" : {
        "amount" : ["1e2", "1e3", "1e4", "5e4"],
        "ref" : "ref-33step.stat"
    }
}
```

### Plotting
You can decide coloumns should be plotted just by adding them to the array called `plotting_cols`. 
Example
```python
plotting_cols = [
    "rms_x",
    "rms_y",
    "rms_s",
    "max_x",
    "max_y",
    "max_s",
    "xpx",
    "ypy",
    "zpz",
    "Dx",
    "DDx",
    "Dy",
    "DDy"
]
```
Note that the coloumn should exist for both opal and opalx stat or else the program will crash

### Template
The input files are being generated using the template.py file. If you want to use different elements change the `get_opalx_string(amount, steps)` function accordingly.

### Cleaning
To remove all stats and plots just run the `clean.sh` bash script. 