This subdirectory contains the reference genome FASTA file
as well the results obtained by FastaToCpG.m. The results
are organized with a subdirectory for each species (Human,
Mouse, etc.), with each species folder including a MATLAB
MAT file CpGlocationChr#.mat for each chromosome that contains
the following information:
o location of CpG sites
o CpG density for each CpG site
o distance between neighboring CpG sites
o location of the last CpG site in the chormosome
o length of chromosome (in # of base pairs)

