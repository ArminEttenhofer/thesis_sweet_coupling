#! /bin/bash


for i in *.bmp; do
	echo "Processing '$i'"
	convert "$i" "${i%.*}.png" #&& rm "$i"
done
