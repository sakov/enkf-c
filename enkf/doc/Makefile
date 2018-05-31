TEXFILES = $(wildcard *.tex)
MAINFILE = enkf-userguide.tex
VIEWER = evince

VIEWER_VERSION := $(shell $(VIEWER) --version 2>/dev/null)
SHELL := /bin/bash -O extglob
TARGETS = $(patsubst %.tex,%.pdf,$(TEXFILES))
MAINTARGET = $(patsubst %.tex,%.pdf,$(MAINFILE))

RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined|Citation(s) may have changed"

all: viewer all-recursive $(TARGETS)
	$(VIEWER) $(MAINTARGET) >/dev/null 2>&1

all-recursive:
	@for dir in $(wildcard *); do if [ -d $$dir ] && [ -f $$dir/Makefile ]; then cd $$dir; $(MAKE) all; cd ..; fi; done

clean: cleanthis
	@for dir in $(wildcard *); do if [ -d $$dir ] && [ -f $$dir/Makefile ]; then cd $$dir; $(MAKE) cleanthis; cd ..; fi; done

cleanthis:
	@rm -f *.aux *.log *.bbl *.blg *.brf *.cb *.ind *.idx *.ilg  \
          *.inx *.ps *.dvi *.toc *.out *.lot *~ *.lof *.ttt *.fff *.spl

cleanall: clean cleanpdf
	@for dir in $(wildcard *); do if [ -d $$dir ]; then cd $$dir; for file in *.tex; do rm -f $${file%.tex}.pdf; done; cd ..; fi; done

cleanpdf:
	@rm -f *.pdf

%.pdf: %.tex
	pdflatex $<
	egrep -c $(RERUNBIB) $*.log && (bibtex $*;pdflatex $<); true
	egrep -c $(RERUNBIB) $*.log && (bibtex $*;pdflatex $<); true
	egrep $(RERUN) $*.log && (pdflatex $<) ; true
	egrep $(RERUN) $*.log && (pdflatex $<) ; true

viewer:
ifndef VIEWER_VERSION
	@echo "  Could not find PDF viewer \"$(VIEWER)\""
	@echo "  Please set variable VIEWER in Makefile"
	@exit 1
endif
