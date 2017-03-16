#################################################
##BEGIN OF makefile "makefile" - using the gnu-compiler g++
##This file is used by the Make command.
###################################################

FILE=opal_user_guide
PREFIX=$(DOC_PREFIX)

all:	$(FILE).tex
	pdflatex -shell-escape $(FILE).tex
	biber --bblsafechars $(FILE)
	makeindex $(FILE).idx
	pdflatex -shell-escape $(FILE).tex

install:
	cp opal_user_guide.pdf $(PREFIX)

clean:
	rm -f *~ *.lot *.bbl *.aux *.blg *.log *.dvi *.ps *.toc *.idx *.lof *.run.xml *.bcf
	rm -f introduction.pdf introduction.out
	rm -f conventions.pdf conventions.out
	rm -f tutorial.pdf tutorial.out
	rm -f opalt.pdf opalt.out
	rm -f opalcycl.pdf opalcycl.out
	rm -f format.pdf format.out
	rm -f control.pdf control.out
	rm -f elements.pdf elements.out
	rm -f lines.pdf lines.out
	rm -f physics.pdf physics.out
	rm -f distribution.pdf distribution.out
	rm -f fieldsolvers.pdf fieldsolvers.out
	rm -f wakefields.pdf wakefields.out
	rm -f geometry.pdf geometry.out
	rm -f track.pdf track.out
	rm -f fieldemission.pdf fieldemission.out
	rm -f multipact.pdf multipact.out
	rm -f partmatter.pdf partmatter.out

introduction:    introduction.tex
	pdflatex -shell-escape introduction.tex
	biber --bblsafechars introduction
	pdflatex -shell-escape introduction.tex

conventions:	conventions
	pdflatex -shell-escape conventsions.tex
	biber --bblsafechars conventsions
	pdflatex -shell-escape conventsions.tex

tutorial:	tutorial
	pdflatex -shell-escape tutorial.tex
	biber --bblsafechars tutorial
	pdflatex -shell-escape tutorial.tex

opalt:		opalt
	pdflatex -shell-escape opalt.tex
	biber --bblsafechars opalt
	pdflatex -shell-escape opalt.tex

opalcycl:	opalcycl
	pdflatex -shell-escape opalcycl.tex
	biber --bblsafechars opalcycl
	pdflatex -shell-escape opalcycl.tex

format:		format.tex
	pdflatex -shell-escape format.tex
	biber --bblsafechars format
	pdflatex -shell-escape format.tex

control:	control.tex
	pdflatex -shell-escape control.tex
	biber --bblsafechars control
	pdflatex -shell-escape control.tex

elements:	elements.tex
	pdflatex -shell-escape elements.tex
	biber --bblsafechars elements
	pdflatex -shell-escape elements.tex

lines:		lines.tex
	pdflatex -shell-escape lines.tex
	biber --bblsafechars lines
	pdflatex -shell-escape lines.tex

physics:	physics.tex
	pdflatex -shell-escape physics.tex
	biber --bblsafechars physics
	pdflatex -shell-escape physics.tex

distribution:	distribution.tex
	pdflatex -shell-escape distribution.tex
	biber --bblsafechars distribution
	pdflatex -shell-escape distribution.tex

fieldsolvers:	fieldsolvers.tex
	pdflatex -shell-escape fieldsolvers.tex
	biber --bblsafechars fieldsolvers
	pdflatex -shell-escape fieldsolvers.tex

wakefields:	wakefields.tex
	pdflatex -shell-escape wakefields.tex
	biber --bblsafechars wakefields
	pdflatex -shell-escape wakefields.tex

geometry:	geometry.tex
	pdflatex -shell-escape geometry.tex
	biber --bblsafechars geometry
	pdflatex -shell-escape geometry.tex

track:		track.tex
	pdflatex -shell-escape track.tex
	biber --bblsafechars track
	pdflatex -shell-escape track.tex

fieldemission:	fieldemission.tex
	pdflatex -shell-escape fieldemission.tex
	biber --bblsafechars fieldemission
	pdflatex -shell-escape fieldemission.tex

multipact:	multipact.tex
	pdflatex -shell-escape multipact.tex
	biber --bblsafechars multipact
	pdflatex -shell-escape multipact.tex

partmatter:      partmatter.tex
	pdflatex -shell-escape partmatter.tex
	biber --bblsafechars partmatter
	pdflatex -shell-escape partmatter.tex