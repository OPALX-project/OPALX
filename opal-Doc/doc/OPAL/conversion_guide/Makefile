#################################################
##BEGIN OF makefile "makefile" - using the gnu-compiler g++
##This file is used by the Make command.
###################################################

FILE=conversion

all:	$(FILE).tex
	latex $(FILE).tex
	latex $(FILE).tex
#       dvips $(FILE).dvi -o $(FILE).ps -x 1372 -O 2mm,5mm
	dvips $(FILE).dvi -o $(FILE).ps -x 1200 -O 10mm,15mm
clean:
	rm -f *~ *.bbl *.aux *.blg *.log *.dvi *.ps *.toc *.idx *.lof *.pdf

