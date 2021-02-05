#!/bin/bash
 
# navigate to IMG_Dump where all compressed genomes were downloaded
cd ./data/genomes/IMG_Dump/

# Create files to move the .gna and .gaa files into
mkdir nt/
mkdir aa/

# loop through and unpack the tarballs 
for file in *.tar.gz
 do
	 tar -zxf "$file"
done

# move all the .genes.fna files into nt/
mv **/*.genes.fna nt/

# concatenate into a single fasta
cd nt/
cat *.fna >> Merged.fna

# Create file lists for fastas
ls *.genes.fna > fileList.txt

# move all the .genes.faa files into aa/
cd ../
mv **/*.genes.faa aa/

# concatenate into a single fasta
cd aa/
cat *.faa >> Merged.faa

# Create file lists for fastas
ls *.genes.faa > fileList.txt