set datafile separator ","
set title 'titletitletitle'
set xlabel "Time (s)"
set ylabel "Velocity ({/Symbol=m}m/min)"
#set y2label "Links Broken"
#set y2tics 0, 10
set border 11
set xtics nomirror
set ytics nomirror
set term postscript eps color enhanced "Helvetica" 12
#set style line 1 lt 1 lw 1
#set style line 3 lt 3 lw 1
#set style line 4 lt 4 lw 1
#set style line 8 lt 8 lw 1
set output "temp.eps"
set size 0.6,0.6
set pointsize 0.7
set key spacing 1.3
#set key autotitle columnhead
set format y "%3.1f"
#set scale factor 0.00001
#set xrange [0:5]
set yrange [0:1.8]
#set y2range [0:60]
set key off

#plot "velocities.txt" using 1:5:(0.0001) with lines smooth acsplines lt 1 notitle axis x1y1, \
#     "velocities.txt" using 1:5 with points lt 1 pt 1 ti "Velocity" axis x1y1

plot "velocities.txt" using 1:($5*60) with lines lt 1 ti "Velocity" axis x1y1

