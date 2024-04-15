#!/usr/bin/env bash

basename=$1
start=$2
step=$3
end=$4
data=$5

gnuplot -e "filename='$basename.gif'" -e "start=$start" -e "end=$end" -e "step=$step" -e "data='$data'" plot.p
##convert $basename.gif -trim +repage $basename.trim.gif
##convert $basename.trim.gif -fuzz 30% -transparent black $basename.trans.gif
##convert forest_map.png -resize 555x394 forest_map_scaledown.png
##convert $basename.trans.gif -coalesce null: forest_map_scaledown.png -compose overlay -layers omposite -layers optimize $basename.overlay.gif
##convert $basename.overlay.gif -coalesce img/$basename.overlay.%03d.png

