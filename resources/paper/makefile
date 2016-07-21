all: nusquids.pdf

TEX_FILES:=nusquids.tex body_track.tex cross_section_tau_decay.tex nusquids_class.tex nusquidsatm_class.tex

nusquids.pdf : $(TEX_FILES) nusquids.bib
	pdflatex nusquids
	bibtex nusquids
	pdflatex nusquids
	pdflatex nusquids

fast : 
	pdflatex nusquids

clean : 
	rm -f nusquids.pdf
	rm -f nusquids.aux nusquids.bbl nusquids.blg nusquids.log nusquids.out nusquids.spl nusquids.toc
	rm -f arxiv.tar.bz2

arxiv : nusquids.pdf
	mkdir -p arxiv
	cp -pr nusquids.tex nusquids.bbl fig arxiv/
	tar cjf arxiv.tar.bz2 arxiv
	rm -rf arxiv

.PHONY : all fast clean arxiv
