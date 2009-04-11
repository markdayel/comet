set datafile separator ","



set term postscript eps color enhanced "Helvetica" 10
set output "temp.eps"

set xzeroaxis lt 3 linecolor rgbcolor "#000000" lw 1
set style line 1 lt 1 linecolor rgbcolor "#A0A0A0" lw 5
set style line 2 lt 1 linecolor rgbcolor "#00974f" lw 1 pt 1
set style line 3 lt 1 linecolor rgbcolor "red" lw 1 pt 6
set style line 4 lt 1 linecolor rgbcolor "cyan" lw 1 pt 2
set style line 5 lt 1 linecolor rgbcolor "orange" lw 1 pt 7

set border 11
set xtics nomirror
set ytics nomirror 
set size 0.35,0.30
set pointsize 0.7

set key spacing 1.0
#set key center right

set origin 1.0,0

#set yrange [0:6]
set xrange [0:5]
set y2range [0:60]
set y2tics 0, 20
set ytics -6, 1

set xlabel "Distance from surface ({/Symbol m}m)" 0.0,0.3
set ylabel "Tensile Force (a.u.)" 1.5,0.0
set y2label "Link Breaks" -0.5,0.0
set key samplen 2  # length of the sample line in the key


# for plotting the radial and tangential forces, summing the link and the repulsive forces 
set yrange [-6:2.5]

set label "titletitletitle" at -12.55,-4.1 rotate left
set label "Tracks" at -11.1,-7
set label "Link Breaks" at -7.5,-7
set label "Circumferential Tension" at -4.6,-7
set label "(i)" at -12.22,2.4 textcolor rgbcolor "white"
set label "(ii)" at -8.5,2.4 textcolor rgbcolor "white"
set label "(iii)" at -4.76,2.4 textcolor rgbcolor "white"
set label "(iv)" at -0.95,2.4 textcolor rgbcolor "black"
set key 5, -1.5

# for plotting the forces separately
# set yrange [0:6]
# set label "titletitletitle" at -12.55,1.5 rotate left
# set label "Tracks" at -11.1,-0.6
# set label "Link Breaks" at -7.5,-0.6 
# set label "Circumferential Tension" at -4.6,-0.6
# set label "(i)" at -12.22,6 textcolor rgbcolor "white"
# set label "(ii)" at -8.5,6 textcolor rgbcolor "white"
# set label "(iii)" at -4.76,6 textcolor rgbcolor "white"
# set label "(iv)" at -0.95,6 textcolor rgbcolor "black"
# set key 5, 5.7


#set nokey
#set key autotitle columnhead
#set format y "%g"
#set scale factor 0.00001


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
# plot "temp.txt" using 1:7 with impulses ls 1 ti "Links Broken" axis x1y2, \
#      "temp.txt" using 1:($3*0.000001):(150.0) with lines smooth acsplines ls 2 notitle axis x1y1, \
#      "temp.txt" using 1:($4*0.000001):(150.0) with lines smooth acsplines ls 3 notitle axis x1y1, \
#      "temp.txt" using 1:($5*0.000001):(150.0) with lines smooth acsplines ls 4 notitle axis x1y1, \
#      "temp.txt" using 1:($6*0.000001):(150.0) with lines smooth acsplines ls 5 notitle axis x1y1, \
#      "temp.txt" using 1:($3*0.000001) with points ls 2  ti "F_{rep}(radial)" axis x1y1, \
#      "temp.txt" using 1:($4*0.000001) with points ls 3  ti "F_{rep}(circ)" axis x1y1, \
#      "temp.txt" using 1:($5*0.000001) with points ls 4  ti "F_{link}(radial)" axis x1y1, \
#      "temp.txt" using 1:($6*0.000001) with points ls 5  ti "F_{link}(circ)" axis x1y1

# plotting the radial and tangential forces, summing the link and the repulsive forces 
	plot "temp.txt" using 1:7 with impulses ls 1 ti "Links Broken" axis x1y2, \
	     "temp.txt" using 1:(($5-$3)*0.000001):(150.0) with lines smooth acsplines ls 2 notitle axis x1y1, \
	     "temp.txt" using 1:(($6-$4)*0.000001):(150.0) with lines smooth acsplines ls 3 notitle axis x1y1, \
	     "temp.txt" using 1:(($5-$3)*0.000001) with points ls 2  ti "F_{radial}" axis x1y1, \
	     "temp.txt" using 1:(($6-$4)*0.000001) with points ls 3  ti "F_{circ}" axis x1y1 \

# for plotting the forces separately
	# plot "temp.txt" using 1:7 with impulses ls 1  ti "Links Broken" axis x1y2, \
	#      "temp.txt" using 1:($3*0.000001):(150.0) with lines smooth acsplines ls 2 notitle axis x1y1, \
	#      "temp.txt" using 1:($4*0.000001):(150.0) with lines smooth acsplines ls 3 notitle axis x1y1, \
	#      "temp.txt" using 1:($5*0.000001):(150.0) with lines smooth acsplines ls 4 notitle axis x1y1, \
	#      "temp.txt" using 1:($6*0.000001):(150.0) with lines smooth acsplines ls 5 notitle axis x1y1, \
	#      "temp.txt" using 1:($3*0.000001) with points ls 2 ti "F_{rep}(radial)" axis x1y1, \
	#      "temp.txt" using 1:($4*0.000001) with points ls 3 ti "F_{rep}(circ)" axis x1y1, \
	#      "temp.txt" using 1:($5*0.000001) with points ls 4 ti "F_{link}(radial)" axis x1y1, \
	#      "temp.txt" using 1:($6*0.000001) with points ls 5 ti "F_{link}(circ)" axis x1y1

#plot "temp.txt" using 1:7 with impulses lt 7 lw 5 ti "Link Breaks" axis x1y2, \
#	     "temp.txt" using 1:($5*0.000001):(150.0) with lines smooth acsplines lt 4 notitle axis x1y1, \
#	     "temp.txt" using 1:($6*0.000001):(150.0) with lines smooth acsplines lt 3 notitle axis x1y1, \
#	     "temp.txt" using 1:($5*0.000001) with points lt 4 pt 2 ti "F_{radial}" axis x1y1, \
#	     "temp.txt" using 1:($6*0.000001) with points lt 3 pt 6 ti "F_{circ}" axis x1y1
