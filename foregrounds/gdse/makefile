target: grf visfitsgrid v2psfitsc avgmgfile binv2psfitscm varv2psfits3dc point

grf:  grf.c fitsprog.c beam.c
	gcc -g -o grf fitsprog.c beam.c grf.c -lfftw3 -lcfitsio -lm -lgsl -lgslcblas
	rm -rf *~
visfitsgrid: fitsprog.c visfitsgrid.c read_fits_func.c 
	gcc -g -o visfitsgrid fitsprog.c visfitsgrid.c read_fits_func.c -lcfitsio -lfftw3 -lm
	rm -rf *~

v2psfitsc: v2psfitsc.c read_fits_func.c
	gcc -o v2psfitsc read_fits_func.c v2psfitsc.c -lcfitsio -lnrcp -lfftw3 -lm
	rm -rf *~

v2psfitsm: v2psfitsm.c read_fits_func.c
	gcc -o v2psfitsm read_fits_func.c v2psfitsm.c -lcfitsio -lnrcp -lfftw3 -lm
	rm -rf *~

binv2psfitscm: binv2psfitscm.c
	gcc -o binv2psfitscm binv2psfitscm.c -lcfitsio -lm
	rm -rf *~

varv2psfits3dc: varv2psfits3dc.c
	gcc -g -o varv2psfits3dc varv2psfits3dc.c -lcfitsio -lm
	rm -rf *~

avgmgfile: avgmgfile.c
	gcc -o avgmgfile avgmgfile.c -lcfitsio -lnrcp -lfftw3 -lm
	rm -rf *~

clean:
	rm -rf grf visfitsgrid v2psfitsc avgmgfile binv2psfitscm varv2psfits3dc
