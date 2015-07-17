# note - this makefile is only for the manuscript and associated files; analyses need to be made by hand

all: manuscript si

bibliography: tex/2state.bib
manuscript: ms_2_state.pdf
figs: img/response_curves.png img/posterior_maps.png
si: si_2_state.pdf

tex/2state.bib: /Users/mtalluto/Dropbox/work/papers/master_db.bib
	bibtally.py tex/manuscript.tex | bibextract.py /Users/mtalluto/Dropbox/work/papers/master_db.bib | bibabbrev.py >tex/2state.bib
	
ms_2_state.pdf: tex/manuscript.tex tex/2state.bib img/posterior_maps.png \
img/response_curves.png tex/pnas.bst
	cd tex; pdflatex manuscript; bibtex manuscript; pdflatex manuscript; pdflatex manuscript
	mv tex/manuscript.pdf ms_2_state.pdf
	
img/posterior_maps.png: scr/fig/fig2_posterior_maps.r
	Rscript scr/fig/fig2_posterior_maps.r
	
img/response_curves.png: scr/fig/fig1_response_curves.r
	Rscript scr/fig/fig1_response_curves.r
	
si_2_state.pdf: tex/si.tex tex/si_tab_col_models.tex tex/si_tab_ext_models.tex tex/si_tab_species.tex
	cd tex; pdflatex si
	mv tex/si.pdf si_2_state.pdf