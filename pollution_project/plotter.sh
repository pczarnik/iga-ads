#!/bin/bash

mkdir -p results

echo run only in the directory with out_xyz.data files
echo generates pngs and creates animation

python3 rename_files.py

for file in out_*.data; do
    base_name=$(basename "$file" .data)
    
    output_image="results/${base_name}.png"
    
    echo "set term png" > plot_script.gp
    echo "set output '$output_image'" >> plot_script.gp
    echo "plot '$file' with image" >> plot_script.gp
    
    gnuplot plot_script.gp
    
    rm plot_script.gp
done

ffmpeg -framerate 30 -pattern_type glob -i 'results/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4
