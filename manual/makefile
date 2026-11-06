TEX = $(wildcard *.tex)
PDF = $(addsuffix .pdf, $(basename $(TEX)))
INPUTBIB = bibliography.bib

BIBTEXFLAGS = 
PDFLATEXFLAGS = -halt-on-error
ISPELLFLAGS = 

default: pdf

pdf: pdf-pdflatex

pdf-pdflatex: $(addsuffix .pdf, $(basename $(TEX))) $(INPUTBIB)

%.pdf: %.tex
	@# Use bash time commmand rather than any system time binary
	bash -c "set -o pipefail; \
	  time { \
	    pdflatex $(PDFLATEXFLAGS) $< | tee $(basename $<).pdflatex.log; \
	    bibtex $(BIBTEXFLAGS) $(basename $<) | tee -a $(basename $<).pdflatex.log; \
	    pdflatex $(PDFLATEXFLAGS) $< | tee -a $(basename $<).pdflatex.log; \
	    pdflatex $(PDFLATEXFLAGS) $< | tee -a $(basename $<).pdflatex.log; \
	    pdflatex $(PDFLATEXFLAGS) $< | tee -a $(basename $<).pdflatex.log; \
	  } \
	"

clean:
	rm -f $(PDF)
	rm -f $(addsuffix .aux, $(basename $(TEX)))
	rm -f $(addsuffix .bbl, $(basename $(TEX)))
	rm -f $(addsuffix .blg, $(basename $(TEX)))
	rm -f $(addsuffix .log, $(basename $(TEX)))
	rm -f $(addsuffix .pdflatex.log, $(basename $(TEX)))
	rm -f $(addsuffix .toc, $(basename $(TEX)))
	rm -f $(addsuffix .out, $(basename $(TEX)))
