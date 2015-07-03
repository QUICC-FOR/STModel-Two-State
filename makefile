# note - this makefile is only for the manuscript and associated files; analysis and figures
# have to be made by hand

bibliography: tex/2state.bib
manuscript: manuscript.pdf

tex/2state.bib: /Users/mtalluto/Dropbox/work/papers/master_db.bib
	bibtally.py tex/manuscript.tex | bibextract.py /Users/mtalluto/Dropbox/work/papers/master_db.bib >tex/2state.bib
	
manuscript.pdf: tex/manuscript.tex tex/2state.bib tex/pnas.bst
	cd tex; pdflatex manuscript; bibtex manuscript; pdflatex manuscript; pdflatex manuscript
	mv tex/manuscript.pdf manuscript.pdf