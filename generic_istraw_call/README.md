# scripts to convert Affymetrix CEL files into genotype calls
These scripts are generic scripts to call genotypes from a batch of Affymetrix
Axiom CEL files. See the two test scripts for examples of calling genotypes
from the IStraw90 and IStraw35 chips. The examples start by getting a list of
input files from a mysql database, but this can be replaced with any suitable
method of generating the input file list.

Requires Affymetrix Power Tools and SNPolisher to be installed.
See the Affymetrix Power Tools and SNPolisher documentation for further details.
