#!/bin/bash

# Assume we have the file id (or not).

# We want to do this for several files.
# Make sure to test this first on 0001.
for ID in 1
do

# The output file.
OutFile="trial_$ID.tex"

# Get the first part of the file.
# Do this if you want to generate a seperate file.
# head -n 9 ~/Latex/template.tex > $OutFile

# Insert the tables.
# Change to append if you want a seperate file.

# For the true values, if we have synthetic data.
cat "Tbl-true-$ID.tex" | \
sed -e '/center/d' -e '/%/d' > $OutFile

# For the mcmc values.
cat "Tbl-mcmc-$ID.tex" | \
sed -e '/center/d' -e '/%/d' >> $OutFile

# Put in a new line.
echo "" >> $OutFile

# For the seed values.
cat "Tbl-seed-$ID.tex" | \
sed -e '/center/d' -e '/%/d' >> $OutFile

# For the prior values.
cat "Tbl-prior-$ID.tex" | \
sed -e '/center/d' -e '/%/d' >> $OutFile

# Put in a new line.
echo "" >> $OutFile

# For the posterior values.
cat "Tbl-post-$ID.tex" | \
sed -e '/center/d' -e '/%/d' >> $OutFile

# Write the rest of the file
echo "" >> $OutFile
echo "\vspace{12pt}" >> $OutFile
echo "\includegraphics{plots_$ID.eps}" >> $OutFile 
echo "" >> $OutFile

# Include this if you want a seperate file.
# echo "\end{document}" >> $OutFile

# Include this if you want a seperate file.
# latex $OutFile

done

# Remove the Tbl-*-$ID.tex files.
# rm "Tbl-*-$ID"

# Now put it all together and make a pdf.
# Again, do this if you want a seperate file.
# dviconcat -o all.dvi trial_*.dvi
# dvipdf all.dvi all.pdf

# Clean things up.
# Again, do this is you want a seperate file.
# rm *.log
# rm *.aux
# rm trial_*.dvi