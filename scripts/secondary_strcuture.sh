#!/bin/bash
f=$(basename $1 .xtc )
echo 1 | gmx do_dssp -f $1 -s $2 -sc ${f}_scount.xvg -a ${f}_area.xpm -o ${f}_ss.xpm -ta ${f}_totarea.xvg -aa ${f}_averarea.xvg
echo 1 | gmx rmsdist -f $1 -s $2 -o ${f}_distrmsd.xvg -rms ${f}_rmsdist.xpm -scl ${f}_rmsscale.xpm -mean ${f}_rmsmean.xpm
