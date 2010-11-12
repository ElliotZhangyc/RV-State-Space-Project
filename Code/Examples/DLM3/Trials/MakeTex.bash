#!/bin/bash

# Assume we have the file id (or not).

# We want to do this for several files.
# Make sure to test this first on 0001.
for ID in $@
do

# The output file.
OutFile="trial_$ID.tex"

# Get the prior data.
prior=$(head "Tbl-prior-$ID.txt")
# Get the posterior data.
post=$(head "Tbl-post-$ID.txt")
# Get the seed data.
seed=$(head "Tbl-seed-$ID.txt")
# Get the mcmc data.
mcmc=$(head "Tbl-mcmc-$ID.txt")
# Get the true values.
true=$(head "Tbl-true-$ID.txt")
# Get the synthetic data.
synth=$(head "Tbl-synth-$ID.txt")

# Insert the data.
cat "table_template.tex" | \
sed -e "s/%prior/$prior/" \
    -e "s/%post/$post/" \
    -e "s/%seed/$seed/" \
    -e "s/%mcmc/$mcmc/" \
    -e "s/%true/$true/" \
    -e "s/%synth/$synth/" \
    > $OutFile

# Write the rest of the file
echo -e "\n" >> $OutFile
echo "\vspace{12pt}" >> $OutFile
echo "\includegraphics{plots_$ID.eps}" >> $OutFile 
echo "" >> $OutFile

echo "Generated file <$OutFile>."

done

# Remove the Tbl-*-$ID.tex files.
# rm "Tbl-*-$ID"