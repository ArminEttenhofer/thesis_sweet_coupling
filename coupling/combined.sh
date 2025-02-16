#! /usr/bin/env bash

precice=""

while getopts p: flag
do
    case "${flag}" in
        p) precice="-p ${OPTARG}";;
        *) ;;
    esac
done


./run.sh ${precice}
./plot.sh
