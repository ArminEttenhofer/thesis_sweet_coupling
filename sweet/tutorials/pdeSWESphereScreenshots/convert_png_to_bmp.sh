#! /bin/bash


for i in *.png; do
	echo "Processing '$i'"
	convert "$i" "${i%.*}.bmp" #&& rm "$i"
done
