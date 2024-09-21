#!/bin/bash
#This script downloads and prepares cage seq data

###
###
###Download bigWigToBedGraph
wget -O  ../../software/ucscTools/bigWigToBedGraph https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod +x ../../software/ucscTools/bigWigToBedGraph

echo "---"
echo "---"
echo "---"
echo "bigWigToBedGraph downloaded"
echo "---"
echo "---"
echo "---"

###
###
###Download plus and minus strand and make bedGraph files
wget -O ../../data/cage/signals/pancreas-Adult.plus.bw https://epd.expasy.org/epd/ucsc/epdHub/hg38/bigWIG/fantom5/pancreas-Adult.plus.bw
wget -O ../../data/cage/signals/pancreas-Adult.minus.bw https://epd.expasy.org/epd/ucsc/epdHub/hg38/bigWIG/fantom5/pancreas-Adult.minus.bw
../../software/ucscTools/bigWigToBedGraph ../../data/cage/signals/pancreas-Adult.plus.bw ../../data/cage/signals/pancreas-Adult.plus.bedGraph
../../software/ucscTools/bigWigToBedGraph ../../data/cage/signals/pancreas-Adult.minus.bw ../../data/cage/signals/pancreas-Adult.minus.bedGraph

echo "---"
echo "---"
echo "---"
echo "CAGE-Seq files downloaded"
echo "---"
echo "---"
echo "---"