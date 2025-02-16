#! /usr/bin/env bash

# Remove old stuff
rm -rf plots

# Run script generation
mkdir plots
cd plots

mkdir h
mkdir u
mkdir v
mkdir mag

cd h
python3 ../../pp_plot.py ../../job*/output_prog_h*
cd ..
ffmpeg -f image2 -r 8 -pattern_type glob -i 'h/output_prog_h_*' -f lavfi -i anullsrc  -c:v libx264 -c:a aac -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -shortest anim_h.mp4

