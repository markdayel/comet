set datafile separator","

set term postscript eps color enhanced "Helvetica" 14
set output "temp.eps"
set size 1.0,.8
set multiplot

#set title 'titletitletitle'

#set border 3
set border 11


#set xtics 10
#set mxtics 2

#set ytics 30
#set mytics 3

set xtics nomirror
set ytics nomirror

#set noxtics
#set noytics


#set style line 1 lt 1 lw 1
#set style line 3 lt 3 lw 1
#set style line 4 lt 4 lw 1
#set style line 8 lt 8 lw 1
set size 0.48,.6
#set pointsize 0.7
#set key spacing 1.3
#set key autotitle columnhead
#set format y "%3.1f"
#set scale factor 0.00001
#set xrange [0:9600]

#set key off


#set grid

set origin 0.5,0
set xlabel "Time (Frame#)"
set ylabel "Stretch Factor"
set y2label "Distance Travelled ({/Symbol m}m)"
set y2tics 0, 5
set xtics 50
set xrange [75:230]
set yrange [0.4:3.2]
set y2range [0:16]

set key off

#plot "velocities.txt" using 1:5:(0.0001) with lines smooth acsplines lt 1 notitle axis x1y1, \
#     "velocities.txt" using 1:5 with points lt 1 pt 1 ti "Velocity" axis x1y1

#plot "dist_mean_sd.txt" using 1:2 with points lt 1 lw 1.5 ti "Angle" axis x1y1
#plot "dist_mean_sd.txt" using 1:2:3 with yerrorbars lt 1 lw 1.5 ti "Angle" axis x1y1

plot "SquashFigData2.csv" using 1:7 with lines lt 6 lc 7 lw 2 ti "Displacement" axis x1y2, \
     "SquashFigData2.csv" using 1:5 with lines lt 1 lc 1 lw 2 ti "Crack" axis x1y1, \
     "SquashFigData2.csv" using 1:3 with lines lt 2 lc 2 lw 2 ti "Shell(inner)" axis x1y1, \
     "SquashFigData2.csv" using 1:4 with lines lt 3 lc 3 lw 2 ti "Shell(outer)" axis x1y1, \
     "SquashFigData2.csv" using 1:6 with lines lt 5 lc 4 lw 2 ti "Radial" axis x1y1, \
     "SquashFigData2.csv" every 15 using 1:7 with points lt 6 lc 7 lw 2 ti "" axis x1y2, \
	  "SquashFigData2.csv" every 15 using 1:5 with points lt 1 lc 1 lw 2 ti "" axis x1y1, \
	  "SquashFigData2.csv" every 15 using 1:3 with points lt 2 lc 2 lw 2 ti "" axis x1y1, \
	  "SquashFigData2.csv" every 15 using 1:4 with points lt 3 lc 3 lw 2 ti "" axis x1y1, \
	  "SquashFigData2.csv" every 15 using 1:6 with points lt 5 lc 4 lw 2 ti "" axis x1y1
     

#     "distsummary1.csv" using 1:5 with lines lt 4 lc 4 lw 2 ti "Outer Crack" axis x1y1, \

set size 0.48,.6

set origin 0,0
set xlabel "Time (seconds)"
set xrange [75:140]
set yrange [0.6:2]
set y2range [0:4]
set y2tics 0, 1
set xtics 25
set key box lw 0.2 150,2.5 reverse top left Left

plot "SquashFigData2.csv" using 8:9:(0.02) with linespoints smooth acsplines lt 6 lc 7 lw 2 ti "Displacement" axis x1y2, \
     "SquashFigData2.csv" using 8:14:(0.03) with linespoints smooth acsplines lt 1 lc 1 lw 2 ti "Crack" axis x1y1, \
     "SquashFigData2.csv" using 8:10:(0.02) with linespoints smooth acsplines lt 2 lc 2 lw 2 ti "Shell(inner)" axis x1y1, \
     "SquashFigData2.csv" using 8:11:(0.02) with linespoints smooth acsplines lt 3 lc 3 lw 2 ti "Shell(outer)" axis x1y1, \
     "SquashFigData2.csv" using 8:12:(0.02) with linespoints smooth acsplines lt 5 lc 4 lw 2 ti "Radial" axis x1y1, \
	 "SquashFigData2.csv" using 8:9 with points lt 6 lc 7 lw 2 ti "" axis x1y2, \
	     "SquashFigData2.csv" using 8:14 with points lt 1 lc 1 lw 2 ti "" axis x1y1, \
	     "SquashFigData2.csv" using 8:10 with points lt 2 lc 2 lw 2 ti "" axis x1y1, \
	     "SquashFigData2.csv" using 8:11 with points lt 3 lc 3 lw 2 ti "" axis x1y1, \
	     "SquashFigData2.csv" using 8:12 with points lt 5 lc 4 lw 2 ti "" axis x1y1

unset multiplot


