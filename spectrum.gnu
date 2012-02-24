set style line 1 lc rgb '#FF0000' lt 1 lw 3 pt 7 ps 1.0
set terminal postscript eps enhanced
set output "spectrum.eps"
set grid
set ytics 1E-8,1E-1,1
set format y "10^{%L}"
set format x "10^{%L}"
unset key
unset mytics
set logscale xy

plot [1:100][1E-8:1] 'spectrum.dat' using 1:2 with linespoints ls 1,\
					 'spectrumSpektralCode.txt'
