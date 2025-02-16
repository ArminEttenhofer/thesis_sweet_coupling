#! /bin/bash

set -e

#for i in ./include/sweet ./programs ./tests ./tutorials; do
for i in ./include/sweet ./programs ./tests; do
	mule.fix_hpp_defines `find $i -name "*.hpp"`
done


