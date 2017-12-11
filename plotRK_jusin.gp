reset

del = gprintf("%.2f",del_);
keisu_omega_a = gprintf("%.2f",keisu_omega_a_);
a = gprintf("%.2f",a_);


#filename = "h10"
datfile = "RK.dat"

set term postscript eps enhanced color "" 25
set output filename

set title "del=".del.", {/Symbol w}_a=".keisu_omega_a."{/Symbol w}_0, a=".a  #font ",80"


set key left top
set xlabel "t" offset 0.0,-1.0 #font ",60"
set ylabel "x(t)" offset -0.0,0.0

set yrange [-pi:pi]


set xtics #font ",60"
set ytics #0.04 #font ",60"

#set key font ",60"

plot datfile \
u 1:3 w l lw 3 lc "red" title "x(t)"  ,\
"" u 1:3 w p notitle,\
"" u 1:4 w l lc "blue" title "x(0)cos({/Symbol w} t)",\
"" u 1:5 w l lc "green" title  "l-acos(w_at - delta)"


unset output
set term x11
