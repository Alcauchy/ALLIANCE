#!/bin/bash

# generate latex files

doxygen alliancedoc

# go to latex source directory and generate a pdf

make -C ./latex
