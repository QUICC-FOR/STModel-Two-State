# note - this makefile is only for the figures; analyses need to be made by hand

figs: img/response_curves.png img/posterior_maps.png
si: si_2_state.pdf

img/posterior_maps.png: scr/fig/fig2_posterior_maps.r
	Rscript scr/fig/fig2_posterior_maps.r
	
img/response_curves.png: scr/fig/fig1_response_curves.r
	Rscript scr/fig/fig1_response_curves.r
	
