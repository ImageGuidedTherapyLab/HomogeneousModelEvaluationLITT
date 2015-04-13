# Makefile for generating presentation $(TARGET).pdf from $(TARGET).tex 
PAPER = gpuPlanningAug13
PRESENTATION = WorkFlow

.PHONY: $(PAPER).pdf 
all: $(PAPER).pdf

# build pdf file first to ensure all aux files are available
# http://latex2rtf.sourceforge.net/
#   M6 converts equations to bitmaps
$(PAPER).rtf: $(PAPER).pdf $(PAPER).bbl $(PAPER).aux
	latex2rtf -M6 -D 600 $(PAPER).tex

$(PRESENTATION).pdf: $(PRESENTATION).tex 
	pdflatex $(PRESENTATION).tex

$(PAPER).pdf: $(PAPER).tex $(PAPER).bbl
	pdflatex $(PAPER).tex

$(PAPER).bbl: $(PAPER).bib 
	pdflatex $(PAPER).tex
	bibtex $(PAPER)
	pdflatex $(PAPER).tex

# build pdf file first to ensure all aux files are available
$(PAPER).dvi: $(PAPER).pdf 
	latex $(PAPER).tex

cover: cover_letter.tex
	pdflatex cover_letter.tex 

clean:
	rm $(PAPER).aux  $(PAPER).dvi $(PAPER).pdf $(PAPER).bbl  $(PAPER).log  $(PAPER).blg  $(PAPER).out

# Environment used for build
# goten$ lsb_release -a
# No LSB modules are available.
# Distributor ID: Ubuntu
# Description:    Ubuntu 12.04.2 LTS
# Release:        12.04
# Codename:       precise
# goten$ dpkg --get-selections | grep tex >> makefile
# gettext						install
# gettext-base					install
# latex-beamer					install
# latex-xcolor					install
# libdjvulibre-text				install
# libexttextcat-data				install
# libexttextcat0					install
# libgettextpo0					install
# liblocale-gettext-perl				install
# libtext-charwidth-perl				install
# libtext-iconv-perl				install
# libtext-wrapi18n-perl				install
# luatex						install
# plymouth-theme-ubuntu-text			install
# preview-latex-style				install
# tex-common					install
# tex4ht						install
# tex4ht-common					install
# texlive-base					install
# texlive-binaries				install
# texlive-common					install
# texlive-doc-base				install
# texlive-extra-utils				install
# texlive-font-utils				install
# texlive-fonts-recommended			install
# texlive-fonts-recommended-doc			install
# texlive-generic-recommended			install
# texlive-latex-base				install
# texlive-latex-base-doc				install
# texlive-latex-extra				install
# texlive-latex-extra-doc				install
# texlive-latex-recommended			install
# texlive-latex-recommended-doc			install
# texlive-luatex					install
# texlive-pictures				install
# texlive-pictures-doc				install
# texlive-pstricks				install
# texlive-pstricks-doc				install
# texlive-science					install
# texlive-science-doc				install
