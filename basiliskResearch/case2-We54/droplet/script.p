unset log
unset label
set xtic auto
set ytic auto
set title "Interface plotting"
set xr [-9.9e-05:0.0006]
set yr [0:0.0005]
plot "interface-0" u 1:2,"interface-7000" u 1:2,"interface-14000" u 1:2,"interface-21000" u 1:2,"interface-28000" u 1:2,"interface-35000" u 1:2