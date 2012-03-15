beginfig(1);
defaultscale:=0.8;
defaultfont:="ptmr8r";
a=0.5cm;
b=0.3cm;
pickup pencircle scaled 3pt;
for i=0 upto 10:
  for j=1 upto 9:
    drawdot(i*a + b,j*b);
  endfor
endfor
pickup pencircle scaled 1pt;
drawarrow (0,0)--(0,9b);
drawarrow (0,0)--(10a+b,0b);
pickup pencircle scaled 0.5pt;
draw (b,b)--(10a+b,b)--(b,2b)--(10a+b,2b)--(b,3b)--(10a+b,3b)--(b,4b)--(10a+b,4b)--(b,5b)--(10a+b,5b)--(b,6b)--(10a+b,6b)--(b,7b)--(10a+b,7b)--(b,8b)--(10a+b,8b)--(b,9b)--(10a+b,9b);
label.bot(btex z (primary direction) etex scaled 0.7, (7a+b,0));
label.lft(btex r (secondary direction) etex rotated 90 scaled 0.7, (0,6b));
endfig;

beginfig(2);
defaultscale:=1.4;
defaultfont:="ptmr8r";
pickup pencircle scaled 3pt;
for i=0 upto 10:
  for j=1 upto 9:
    drawdot(i*a + b,j*b);
  endfor
endfor
pickup pencircle scaled 1pt;
drawarrow (0,0)--(0,9b);
drawarrow (0,0)--(10a+b,0b);
pickup pencircle scaled 0.5pt;
draw (b,b)--(b,9b)--(a+b,b)--(a+b,9b)--(2a+b,b)--(2a+b,9b)--(3a+b,b)--(3a+b,9b)--(4a+b,b)--(4a+b,9b)--(5a+b,b)--(5a+b,9b)--(6a+b,b)--(6a+b,9b)--(7a+b,b)--(7a+b,9b)--(8a+b,b)--(8a+b,9b)--(9a+b,b)--(9a+b,9b)--(10a+b,b)--(10a+b,9b);
label.bot(btex z (secondary direction) etex scaled 0.7, (7a+b,0));
label.lft(btex r (primary direction) etex rotated 90 scaled 0.7, (0,6b));
endfig; 
end
