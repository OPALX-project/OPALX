#!/bin/bash
for f do
    sed "s/writeError(/writeError(errModeAll,"/g $f > tmp.txt
    mv tmp.txt $f
done
