In this directory, scripts for analyses with reference sequences of C. elegans ce11 and E. coli REL606 are located.

# how to use

1. create a directory for a sample and a parameter set
2. move to the directory and create a setting file for the location of subreads data for SMRT Link (.subreads.xml) with command `dataset create`.
3. call alignment and IPD (inter-pulse duration) summarization pipeline with `run_pbsmrtpipe.sh`
4. call extra scripts with `make -f ../samples.makefile`
