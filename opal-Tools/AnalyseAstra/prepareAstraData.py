#!/usr/bin/python

import os, re, sys, string, math

if not len(sys.argv) >= 2 or not os.path.exists(sys.argv[1]):
    print "usage: " + sys.argv[0] + " inputfile.in"
    sys.exit()

inputfn = sys.argv[1]
basefn = inputfn.split('.')[0]
inp = open(inputfn, "r")

runstr = '000'
for lstr in inp:
    lstr = lstr.lstrip()
    run = re.search("^RUN=(\d*)", lstr)
    if run:
        runstr = "%03i" % int(run.group(1))

pwd = os.getenv("PWD")
psdumps = [] ## file numbers (\d{4}) of phase space dumps between 2.95 m and 3.1m
dumppoints = [] ## array of position of ref particle in phase space dumps
refpoints = [] ## postion of ref particle in test.track.004
mass = 510998.910 ## electron rest mass

if not os.path.exists("astra_data"):
    os.mkdir("astra_data")
elif os.path.exists("astra_data/astra.out"):
    os.remove("astra_data/astra.out")

## add z and pz of the ref particle to the phase space dumps
## and convert the momentum to \beta\gamma
## input file(s): basefn.\d{4}.run
## output file(s): astra_data/basefn.\d{4}.run.out
for file in os.listdir(pwd):
    psdump = re.search("^" + basefn + "\.(\d{4})\." + runstr + "$", file)
    if psdump:
        inp = open(file, "r")
        lines = inp.readlines()
        inp.close()
        
        z = float(lines[0].split()[2]);
        if z >= 2.95 and z <= 3.1:
            psdumps.append(psdump.group(1))
            dumppoints.append(z)
            pz = float(lines[0].split()[5])
            out = open("astra_data/" + file + ".out", "w")
            line = lines[0].split()
            out.write(line[0] + "\t" + \
                      line[1] + "\t" + \
                      line[2] + "\t" + \
                      line[3] + "\t" + \
                      line[4] + "\t" + \
                      str(float(line[5]) / mass) + "\n")

            for lstr in lines[1:]:
                line = lstr.split()
                out.write(line[0] + "\t" + \
                          line[1] + "\t" + \
                          str(z + float(line[2])) + "\t" + \
                          line[3] + "\t" + \
                          line[4] + "\t" + \
                          str((pz + float(line[5])) / mass) + "\n")
            out.close()

if len(psdumps) == 0:
    print "Did you run Astra? Could not find any file '" + basefn + ".\\d{4}." + runstr + "'"
    sys.exit()
dumppoints.sort()
psdumps.sort()

## combine all the phase space dumps into one file
## input file(s): astra_data/basefn.\d{4}.run.out
## output file(s): astra_data/astra.out
out = open("astra_data/astra.out", "w")
for file in psdumps:
    inp = open("astra_data/" + basefn + "." + file + "." + runstr + ".out", "r")
    for line in inp:
        out.write(line)
    inp.close()
out.close()


## read in the longitudinal field of the particles, then integrate them
## to calculate the kinetic energy for every particle at every position.
## The initial kinetic energy for the particles are taken from the file
## with the lowest number, xxxx.  From the energy calc the longitudinal
## momentum. From all the positions in the input file choose for every
## phase space dump those two which are closest to it and interpolate
## linearly the position and energies. Repeate this last step also for
## the input file
## input file(s): basefn.track.run
##                astra_data/basefn.xxxx.run.out
## output file(s): astra_data/basefn.track.run.out
##                 astra_data/astra.out.rec
inp = open(basefn + ".track." + runstr, "r")
lines = inp.readlines()
inp.close()

for lstr in lines:
    if  int(lstr.split()[0]) == 1:
        refpoints.append(float(lstr.split()[2]))
        
inp = open("astra_data/" + basefn + "." + psdumps[0] + "." + runstr + ".out","r")
energies = []
## read in the initial momenta and calculate the initial kinetic energies
for line in inp:
    momentum = float(line.split()[5])
    energies.append((math.sqrt(1. + momentum**2) - 1.) * mass)
inp.close()

Np = int(lines[0].split()[0]);
F = []  ## will store the positions and field values for every particle separately
G = []  ## will store the kinetic energies for every particle separately
for i in range(0, Np):
    F.append([[],[]])

for lstr in lines:
    line = lstr.split()
    index = int(line[0]) - 1
    F[index][0].append(float(line[2]))
    F[index][1].append(float(line[5]))

## disregard those positions where the field should be zero (in fact in astra
## they're not zero) in front of the TWS
for k in range(0, len(F[0][0])):
    if F[0][0][k] >= 2.95:
        break

for i in range(0, Np):
    G.append([energies[i]])
    for j in range(1,k+1):
        G[i].append(energies[i])
    ## integrate the field. Since we are dealing with electrons the integral has
    ## to be subtracted from the initial energy
    for j in range(k+1,len(F[i][0])):
        G[i].append(G[i][j-1] - (F[i][0][j] - F[i][0][j-1]) * 0.5 * (F[i][1][j] + F[i][1][j-1]))

    ## calculate the momentum from the kinetic energy
    for j in range(0, len(G[0])):
        G[i][j] = math.sqrt((G[i][j] / mass + 1.)**2 - 1)

lines2 = []
## out = open("test.dat", "w")
for j in range(0, len(G[0])):
    for i in range(Np, 0, -1):
        lines2.append([F[i-1][0][j], G[i-1][j]])
##         out.write("%2.12e\t%2.12e\t%i\n" % (F[i-1][0][j], G[i-1][j], i))
## out.close()

out = open("astra_data/" + basefn + ".track." + runstr + ".dat","w")
out2 = open("astra_data/astra.out.rec", "w")
init_index = 0
j = 0
for point in dumppoints:
    diff = math.fabs(point - refpoints[init_index])
    i = init_index + 1
    while i < len(refpoints):
        next_diff = math.fabs(point - refpoints[i])
        if diff >= next_diff:
            diff = next_diff
            i += 1
        else:
            init_index = i
            i -= 1
            diff = math.fabs(point - refpoints[i-1])
            next_diff = math.fabs(point - refpoints[i+1])
            j = i - 1
            if next_diff < diff:
                j = i + 1
            break

    tau = 0.5;
    if refpoints[j] != refpoints[i]:
        tau = (point - refpoints[i]) / (refpoints[j] - refpoints[i])
##     print "%2.5e\t%2.5e\t%2.5e\t%2.5e" % (refpoints[i], point, refpoints[j],tau)
        
    for k in range(0,Np):
        ivalues = lines[Np * i + k].split()
        jvalues = lines[Np * j + k].split()
        out.write("%2.12e\t%2.12e\t%s\t%s\n" % ((1. - tau) * float(ivalues[2]) + tau * float(jvalues[2]), \
                                                (1. - tau) * float(ivalues[5]) + tau * float(jvalues[5]),
                                                ivalues[0], jvalues[0]))

    for k in range(0,Np):
        ivalues = lines2[Np * i + k]
        jvalues = lines2[Np * j + k]
        out2.write("%2.12e\t%2.12e\t%s\t%s\n" % ((1. - tau) * float(ivalues[0]) + tau * float(jvalues[0]), \
                                                 (1. - tau) * float(ivalues[1]) + tau * float(jvalues[1]),
                                                 ivalues[0], jvalues[0]))


out.close()
out2.close()


commands = """set xlabel 's [m]'
set ylabel '{/Symbol b}_z{/Symbol g} [1]'
set y2label 'E_z [MV/m]'
set xrange [2.95:3.1]
set yrange [13.7:17.7]
set y2range [-2e7:2e6]
set ytics nomirror 
set y2tics
set mytics
set my2tics
set grid
plot 'astra_data/astra.out' u 3:6 w p pt 6 lc 1 t 'Astra orig'
repl 'astra_data/astra.out.rec' u 1:2 w p pt 6 lc 3 t 'Astra reconstr'
"""

commands += "repl 'astra_data/" + basefn +  ".track." + runstr + ".dat' u 1:2 w p pt 4 lc 1 axis x1y2 t 'Astra E_z'\n"
commands += """pause -1 'hit return to continue'
set xrange [2.95:3.03]
set yrange [13.7:15.4]
repl
pause -1 'hit return to continue'
"""

out = open("astra_data/gnu_commands.gpl","w")
out.write(commands)
out.close()
