#! /usr/bin/env bash

rm -rf plots
mkdir plots
cd plots

mkdir u
mkdir v
mkdir mag

python3 ../pp_plot_plane_unstablejet.py ../output/output_{prog_u,prog_v}_*

mkdir h
cd h
python3 ../../pp_plot_plane_unstablejet.py ../../output/output_prog_h_pert*
cd ..


ffmpeg -f image2 -r 8 -pattern_type glob -i 'h/prog_h_*' -f lavfi -i anullsrc  -c:v libx264 -c:a aac -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -shortest anim_h.mp4
ffmpeg -f image2 -r 8 -pattern_type glob -i 'u/prog_u_*' -f lavfi -i anullsrc  -c:v libx264 -c:a aac -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -shortest anim_u.mp4
ffmpeg -f image2 -r 8 -pattern_type glob -i 'mag/prog_mag_*' -f lavfi -i anullsrc  -c:v libx264 -c:a aac -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -shortest anim_mag.mp4
ffmpeg -f image2 -r 8 -pattern_type glob -i 'v/prog_v_*' -f lavfi -i anullsrc  -c:v libx264 -c:a aac -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -shortest anim_v.mp4