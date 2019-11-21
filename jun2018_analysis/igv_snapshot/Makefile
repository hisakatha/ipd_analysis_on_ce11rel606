svgs := $(wildcard *.svg)

pdfs := $(svgs:.svg=.pdf)
pngs := $(svgs:.svg=.png)

all: $(pngs)

# Generated PDFs are too big in file size. PNGs seems to be good.
%.pdf: %.svg
	inkscape --without-gui --export-area-drawing --export-pdf=$@ $<

%.png: %.svg
	inkscape --without-gui --export-area-drawing --export-png=$@ -d 180 $<

