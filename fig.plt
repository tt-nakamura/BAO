reset
set terminal postscript eps enhanced 24
set output "fig.eps"
set key left Left reverse
set xrange [9e-3:4e-1]
set yrange [5e-3:5]
set ylabel 'power spectrum of density fluctuation'
set xlabel 'wave number  / h/Mpc'
set logscale xy
set format xy "10^{%T}"
set label "{/Symbol W}_d = 23.3%" at gr 0.72,0.5
set label "{/Symbol W}_b = 4.63%" at gr 0.72,0.43
set label "{/Symbol W}_{/Symbol L} = 72.1%" at gr 0.72,0.36
set label "A = 7.8x10^{-10}" at gr 0.72,0.29
set label "b = 1.6" at gr 0.72,0.22
b = 1.6
A = 7.8e-10*b**2
plot\
	'Anderson_2013_CMASSDR11_power_spectrum_pre_recon.txt' u 1:($2*$1**3/(2*pi*pi)) t 'observed by SDSS3' w p pt 1,\
	'bao.txt' u 1:($2*A) t 'predicted by Big-Bang theory' w l lt 1
