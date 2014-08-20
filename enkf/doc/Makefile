TEXFILES = enkf-userguide.tex
VIEWER = evince

SHELL := /bin/bash -O extglob
TARGETS = $(patsubst %.tex,%.pdf,$(TEXFILES))

RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined|Citation(s) may have changed"

all: all-recursive $(TARGETS)
	$(VIEWER) $(TARGETS) >/dev/null 2>&1

clean: clean-recursive
	rm -f *.aux *.log *.bbl *.blg *.brf *.cb *.ind *.idx *.ilg  \
          *.inx *.ps *.dvi *.toc *.out *.lot *~ *.lof *.ttt *.fff *.spl

all-recursive:
	for dir in $(wildcard *); do if [ -d $$dir ] && [ -f $$dir/Makefile ]; then cd $$dir; $(MAKE) all; cd ..; fi; done

clean-recursive:
	for dir in $(wildcard *); do if [ -d $$dir ] && [ -f $$dir/Makefile ]; then cd $$dir; $(MAKE) clean; cd ..; fi; done

cleanall:
	rm -f !(*.pdf)

%.pdf: %.tex
	pdflatex $<
	egrep -c $(RERUNBIB) $*.log && (bibtex $*;pdflatex $<); true
	egrep -c $(RERUNBIB) $*.log && (bibtex $*;pdflatex $<); true
	egrep $(RERUN) $*.log && (pdflatex $<) ; true
	egrep $(RERUN) $*.log && (pdflatex $<) ; true
