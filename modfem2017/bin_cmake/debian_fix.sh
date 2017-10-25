#!/bin/bash

# Here is fix for debian based systems and MKL (without this, no binary will ever built because of linkage failure).
# Date: 04.11.2015
# Author: kamich@agh.edu.pl
if [ -f /etc/debian_version ]; then
	echo " Altering linker input"
	linker_files="$(find ./ -name 'link.txt')"
	echo " Found linker files: $linker_files. Now updateing "
	old="libmkl_core.a -Wl,--end-group"
	new="libmkl_core.a  -lm -ldl -Wl,--end-group"
	replaced="$(sed -i "s/$old/$new/g" $linker_files)"	  
fi
#end of fix.
