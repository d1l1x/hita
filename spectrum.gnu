set style line 1 lc rgb '#0060ad' lt 1 lw 3 pt 7 ps 1.0
set terminal postscript
set output "spectrum.eps"
unset key
set logscale x
set logscale y

plot 'spectrum.dat' using 1:2 with linespoints ls 1
