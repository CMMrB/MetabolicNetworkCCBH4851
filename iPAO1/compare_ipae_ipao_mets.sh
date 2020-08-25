#!/bin/bash

for i in $(cat ipae.sbml | grep "species id"| cut -d "\"" -f 2-2| grep -v "_e"|cut -d "_" -f 2-2)
do
    echo "##########################################################"
    echo "MET $i"
    echo "iPAE:"
    echo "$(cat ipae.sbml | grep "species id"| grep "$i" | grep -v "_e")"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "iPAO:"
    echo "$(cat iPAO1_mets.csv | grep "$i")"
    echo "---------------------------------------------------------"
done
