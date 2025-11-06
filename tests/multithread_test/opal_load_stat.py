import re
import pandas as pd

# this pattern checks for the coloumns
COL_PATTERN = r"&column[\S\s]*?name=(.*),[\S\s]*?description=\"\d+"
COL_PATTERN_REGEX = re.compile(COL_PATTERN)

def load_dataset(filename: str):
    lines = ""
    with open(filename) as f:
        lines = f.read()
    
    cols = COL_PATTERN_REGEX.findall(lines)
    # find the begining of the data 
    startofdata = lines.split("&data")[-1].split("end")[-1].strip().splitlines()[lines.count("&parameter"):]
    #amount_of_rows = len(startofdata)
    startofdata = "\n".join(startofdata)
    
    from io import StringIO
    return pd.read_csv(StringIO(startofdata), sep="\\s+", header=None, names=cols)
    


if __name__ == "__main__":
    data = load_dataset("opalx/run1-1e7-6.stat")
    print(data)


