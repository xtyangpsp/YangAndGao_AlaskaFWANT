This folder includes the source data and the plotting script used in Figure 3 in the supplementary information.

Formation of the phase delay data files:
1. The files are named by the iteration number. For example, ite_0.05deg_01_measure_result_Alaska.dat_latlon is the data for the 1st iteration using the reference velocity model, for grid size of 0.05 degrees in both latitude and longitude directions.
2. The data file includes 11 columns.
Column 	1: virtual source name
	2: receiver name
	3-4: latitude and longitude (0-360 format) of the virtual source
	5-6: latitude and longitude (0-360 format) of the receiver
	7: phase delay measured between the observed empirical Green's functions and the synthetics
	8: error of the phase delay measurements, computed from the variation of phase delays among all monthly stacks of the empirical Green's functions;
	9: frequency label with index (f1-f8). corresponding period bands in this study: pband=[100, 200; 75,150;50,100;35,75;25,50;15,35;10,25;7.5,15].
	10. Correlation coefficients for the measurement between the observed data and the synthetics;
	11: signal to noise ratio of the observed empirical Green's function

5C.HRLQ II.KDAK 59.492400 221.171000 57.782800 207.417000 0.198000 1.330000 f6 0.740000 5.600000
5C.HRLQ II.KDAK 59.492400 221.171000 57.782800 207.417000 0.566000 1.340000 f7 0.740000 5.800000
5C.HRLQ IU.COLA 59.492400 221.171000 64.873599 212.138000 0.423000 1.820000 f6 0.640000 6.500000
...

Prepared by Xiaotao Yang 
Contact: stcyang@gmail.com

Reference:
Yang, X., H. Gao (in revision). Seismic imaging of slab segmentation and 
correlation with volcano distribution along the Aleutian-Alaska subduction zone, 
Nature Communications