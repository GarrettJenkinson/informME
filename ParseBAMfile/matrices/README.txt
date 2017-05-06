This subdirectory contains methylation matrices, obtained
by processing available BAM files, organized with a
subdirectory for each species (Human, Mouse, etc.), with
each species folder including subdirectories for each
chromosome (chr1, chr2, etc.). Each chromosome subdirectory
includes the MATLAB file BAMfileName.mat for each BAM
file processed, with each file containing the following
information for each genomic region to be subsequently
used in model estimation:
o data matrix with -1,0,1 values for methylation status
o CpG locations broken down by region.


