#################################################
##BEGIN OF makefile "makefile" - using the gnu-compiler g++
##This file is used by the Make command.
###################################################

FILE=opal_user_guide

all:	$(FILE).tex
	pdflatex $(FILE).tex
	makeindex $(FILE).idx
	pdflatex $(FILE).texa
	pdflatex $(FILE).texa
	cp opal_user_guide.pdf /afs/psi.ch/project/amas/www/docs/opal/opal_user_guide-1.1.8.pdf


clean:
	rm -f *~ *.lot *.bbl *.aux *.blg *.log *.dvi *.ps *.toc *.idx *.lof *.pdf