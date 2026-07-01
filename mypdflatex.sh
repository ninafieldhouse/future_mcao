#!/bin/bash
Usage="Usage: mypslatex2pdf filename (without .tex or any extension)"
if [ $# != 1 ] ; then
  echo $Usage
  exit 1
fi
pdflatex $1.tex
pdfcrop $1.pdf
mv $1-crop.pdf $1.pdf
rm $1.aux $1.log
