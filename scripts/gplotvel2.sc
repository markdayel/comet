set datafile separator ","
#set title 'titletitletitle'
set xlabel "Time (s)"
set ylabel "Velocity ({/Symbol=m}m/min)"
#set y2label "Links Broken"
#set y2tics 0, 10
set border 3

set xtics 1000
set mxtics 5

set ytics 1
set mytics 5

set xtics nomirror
set ytics nomirror

#set noxtics
#set noytics
set term postscript eps color enhanced "Helvetica" 8
#set style line 1 lt 1 lw 1
#set style line 3 lt 3 lw 1
#set style line 4 lt 4 lw 1
#set style line 8 lt 8 lw 1
set output "temp.eps"
set size 0.3,0.2
set pointsize 0.7
set key spacing 1.3
#set key autotitle columnhead
set format y "%3.1f"
#set scale factor 0.00001
set xrange [0:5600]
#set xrange [0:11200]
set yrange [0:1.6]
#set y2range [0:60]
set key off


#plot "velocities.txt" using 1:5:(0.0001) with lines smooth acsplines lt 1 notitle axis x1y1, \
#     "velocities.txt" using 1:5 with points lt 1 pt 1 ti "Velocity" axis x1y1

plot "velocities.txt" using 1:($5*60) with lines lt 1 lw 1.5 ti "Velocity" axis x1y1

