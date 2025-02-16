#! /bin/bash


OUTPUTFILE="output_animation_smaller.gif"

rm -f "${OUTPUTFILE}"

# -loop 0: Loop forever
convert \
	-delay 10 \
	-loop 0 \
	-trim +repage	\
	-resize 400x	\
	*0[2-9]?[0].bmp	\
	*[1-9]??[0].bmp	\
	"${OUTPUTFILE}"
