#!/usr/bin/env gnuplot

set terminal gif size 556,395 animate delay 5 loop 1

unset border
unset xtics
unset ytics
unset raxis
unset rtics
unset key
unset colorbox

set lmargin at screen 0
set rmargin at screen 1
set tmargin at screen 0
set bmargin at screen 1

#set cbrange [0:1]

set output filename

set palette rgb 34,35,36

do for [i=start:end:step] {
  #set title 'time '.i
  plot data.'_'.i.'.data'  with image
}

