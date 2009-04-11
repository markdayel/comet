set datafile separator ","
set title 'titletitletitle'
set xlabel "Distance from surface ({/Symbol m}m)"
set ylabel "Force"
set y2label "Links Broken"
set y2tics 0, 10
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
#set format y "%g"
#set scale factor 0.00001
set xrange [0:5]
set yrange [0:40]
set y2range [0:60]
#plot "temp.txt" using 1:7 with impulses lt 7 lw 5 axis x1y2, \
#     "temp.txt" using 1:3 with linespoints lt 1 pt 7 axis x1y1, \
#     "temp.txt" using 1:4 with linespoints lt 2 pt 6 axis x1y1, \
#     "temp.txt" using 1:5 with linespoints lt 3 pt 1 axis x1y1, \
#     "temp.txt" using 1:6 with linespoints lt 4 pt 2 axis x1y1
#plot "temp.txt" using 1:7 with impulses lt 7 lw 5 axis x1y2, \
#     "temp.txt" using 1:3 with lines lt 1 axis x1y1, \
#     "temp.txt" using 1:4 with lines lt 2 axis x1y1, \
#     "temp.txt" using 1:5 with lines lt 3 axis x1y1, \
#     "temp.txt" using 1:6 with lines lt 4 axis x1y1
#plot "temp.txt" using 1:7 with impulses lt 7 lw 5 axis x1y2, \
#     "temp.txt" using 1:3 with lines lt 1 axis x1y1, \
#     "temp.txt" using 1:4 with lines lt 2 axis x1y1, \
#     "temp.txt" using 1:5 with lines lt 3 axis x1y1, \
#     "temp.txt" using 1:6 with lines lt 4 axis x1y1
#plot "temp.txt" using 1:7 with impulses lt 7 lw 5 axis x1y2, \
#     "temp.txt" using 1:3 with linespoints lt 1 pt 7 axis x1y1, \
#     "temp.txt" using 1:4 with linespoints lt 2 pt 6 axis x1y1, \
#     "temp.txt" using 1:5 with linespoints lt 3 pt 1 axis x1y1, \
#     "temp.txt" using 1:6 with linespoints lt 4 pt 2 axis x1y1
plot "temp.txt" using 1:7 with impulses lt 7 lw 5 ti "Links Broken" axis x1y2, \
     "temp.txt" using 1:($8*0.000001):(150.0) with lines smooth acsplines lt 1 notitle axis x1y1, \
     "temp.txt" using 1:($9*0.000001):(150.0) with lines smooth acsplines lt 2 notitle axis x1y1, \
     "temp.txt" using 1:($10*0.000001):(150.0) with lines smooth acsplines lt 3 notitle axis x1y1, \
     "temp.txt" using 1:($11*0.000001):(150.0) with lines smooth acsplines lt 4 notitle axis x1y1, \
     "temp.txt" using 1:($8*0.000001) with points lt 1 pt 1 ti "E_{rep}(radial)" axis x1y1, \
     "temp.txt" using 1:($9*0.000001) with points lt 2 pt 7 ti "E_{rep}(trans)" axis x1y1, \
     "temp.txt" using 1:($10*0.000001) with points lt 3 pt 2 ti "E_{link}(radial)" axis x1y1, \
     "temp.txt" using 1:($11*0.000001) with points lt 4 pt 6 ti "E_{link}(trans)" axis x1y1

