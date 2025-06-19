import os
import shutil
import subprocess

ROOT_FOLDER = "out"
SLURM_FILE_NAME = "run.slurm"
INPUT_FILE = "template.in"
OPALX_EXECUTABLE_FILE = "/data/user/binder_j/opalx-elements/build/src/opalx"

# beatiful colors
def color(r, g, b, background = False):
    return f"\033[38;2;{r};{g};{b}m"

RED = color(255, 0, 0)
MAGENTA = color(245, 66, 242)
BLUE = color(0, 0, 255)
GREEN = color(0, 255, 0)
ORANGE = color(255, 198, 77)
WHITE = "\033[0m"

BOLD="\033[1m"
RESET="\033[0m"

def get_slurm_string(nodes, threads, exe):
    return f"""#!/bin/bash
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes={nodes}               # No. of nodes
#SBATCH --ntasks-per-node=1     # No. of MPI ranks per node. Merlin CPU nodes have 44 cores
#SBATCH --cpus-per-task={threads}      # No. of OMP threads
#SBATCH --time=00:30:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=nomultithread    # Without hyperthreading
#SBATCH --exclusive            # The allocations will be exclusive if turned on (remove extra hashtag to turn on)

#SBATCH --output={nodes}-{threads}.out  # Name of output file
#SBATCH --error={nodes}-{threads}.err    # Name of error file

export OMP_NUM_THREADS={threads}
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


srun --cpus-per-task={threads} {exe} {INPUT_FILE} --info 10
    """

parameters = [
    # nodes, threads
    [ 1    , 1],
    [ 1    , 2],
    [ 1    , 4],
    [ 1    , 8],
    [ 1    , 16],
    [ 1    , 32],
    [ 2    , 32],
    [ 4    , 32],
]

def get_folder_name(n, t):
    return f"n{n}t{t}"

def generate_folder():
    for (n, t) in parameters:
        os.makedirs(os.path.join(ROOT_FOLDER, get_folder_name(n, t)), exist_ok = True)

def gen_slurm_inputfile():
    for (n, t) in parameters:
        file_path = os.path.join(ROOT_FOLDER, get_folder_name(n, t))

        with open(os.path.join(file_path, SLURM_FILE_NAME), "w") as file:
            file.write(get_slurm_string(n, t, OPALX_EXECUTABLE_FILE))

        shutil.copy(INPUT_FILE, file_path)

def run_slurm():
    os.chdir(ROOT_FOLDER)
    for (n, t) in parameters:
        os.chdir(get_folder_name(n, t))
        try:
            print(f"Starting {ORANGE}Slurm script{WHITE} in folder {BLUE}{get_folder_name(n, t)}{WHITE}", end="", flush=True)

            result = subprocess.run(["sbatch", "run.slurm"], capture_output=True, text=True, check=True)

            if result.returncode == 0:
                print(f" --- {GREEN}Done{WHITE}")
            else:
                print(f"{RED}Return Code: {result.returncode + WHITE}")
                print(result.stdout)

        except subprocess.CalledProcessError as e:
            print(f"Command failed with return code {e.returncode}:")
            print(f"Standard Error:\n{e.stderr}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        os.chdir("..")
    os.chdir("..")


if __name__ == "__main__":
    generate_folder()
    gen_slurm_inputfile()
    run_slurm()