#!/bin/bash

excludeList='Classic/Algorithms/CoordinateSystemTrafo.h Classic/Algorithms/CoordinateSystemTrafo.cpp Structure/BoundingBox.h Structure/BoundingBox.cpp AbsBeamline/ElementBase.h AbsBeamline/ElementBase.cpp AbsBeamline/Component.h AbsBeamline/Component.cpp'

for f in `grep -l -R "PartBunch<double, 3>" AbsBeamline`
 do
     if [[ " $excludeList " =~ .*\ ${f:2}\ .* ]]; then
         echo "skipp $f"
     else
         change -f "PartBunch<double, 3>" "PartBunch_t" $f
     fi

 done
exit 0
