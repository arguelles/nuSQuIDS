all: nusquids.pdf

TEX_FILES:=nusquids.tex
#body_track.tex cross_section_tau_decay.tex nusquids_class.tex nusquidsatm_class.tex

nusquids.pdf : $(TEX_FILES) nusquids.bib nusquids.bbl
	pdflatex nusquids
	pdflatex nusquids

nusquids.bbl : nusquids.tex nusquids.bib
	pdflatex nusquids
	bibtex nusquids

nusquids_journal.bbl : nusquids_journal.tex nusquids.tex nusquids.bib
	pdflatex nusquids_journal
	bibtex nusquids_journal

nusquids_journal.pdf : nusquids_journal.tex nusquids.tex nusquids_journal.bbl
	pdflatex nusquids_journal
	pdflatex nusquids_journal

manual.bbl : manual.tex manual.bib
	pdflatex manual
	bibtex manual

manual.pdf : manual.tex manual.bbl nusquids.tex
	pdflatex manual
	pdflatex manual

fast : 
	pdflatex nusquids

clean : 
	rm -f nusquids.pdf
	rm -f nusquids.aux nusquids.bbl nusquids.blg nusquids.log nusquids.out nusquids.spl nusquids.toc
	rm -f nusquids_journal.pdf
	rm -f nusquids_journal.aux nusquids_journal.bbl nusquids_journal.blg nusquids_journal.log nusquids_journal.out nusquids_journal.spl nusquids_journal.toc
	rm -f manual.pdf
	rm -f manual.aux manual.bbl manual.blg manual.log manual.out manual.spl manual.toc
	rm -f arxiv.tar.bz2

arxiv : nusquids.pdf
	mkdir -p arxiv
	cp -pr nusquids.tex nusquids.bbl fig arxiv/
	tar cjf arxiv.tar.bz2 arxiv
	rm -rf arxiv

.PHONY : all fast clean arxiv
