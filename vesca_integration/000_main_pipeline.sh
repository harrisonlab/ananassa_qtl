#!/bin/bash

#
# reconsile the vesca position of markers to the consensus map
#

#select marker blocks which are on the wrong chromosome
#shift them to the end of their correct chromosome
correct_chromosome.R

#select break points, automatically invert and shuffle the blocks as required into optimal order
#including blocks appended to chromosome ends previously
autofix_vesca_posns.R

