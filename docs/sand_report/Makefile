SHELL = /bin/sh

master = main.tex
target = $(master:.tex=)

files = *.tex *.bib 
destbib = references.bib

default: pdf
all: pdf

pdf: $(master:.tex=.pdf)

clean:
	@for t in aux log nav out pdf snm toc dvi bbl blg spl lot lof brf ; do \
	  /bin/rm -f $(target).$$t ; \
	done
	@/bin/rm -f core
$(master:.tex=.pdf): $(files)
	pdflatex $(target) -shell-escape
	@if [ -e $(target).aux ]; then \
	  \
	  bibtex $(target) ;\
	  pdflatex $(target) ;\
	  pdflatex $(target) ;\
	fi 
	pdflatex $(target)
