#!/bin/bash

excludeList=''

for f in `grep -l -R "PartBunch_t" `
 do
     if [[ " $excludeList " =~ .*\ ${f:2}\ .* ]]; then
         echo "skipp $f"
     else
         change -f  "PartBunch_t" "PartBunch_t"  $f
         #change -f "endl\)\;" "endl\;" $f
#         emacs $f

     fi

 done
exit 0
