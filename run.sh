#!/bin/bash
./compile.sh
./execargmc.out > output.txt
python prettyplot.py
#open prettyplot.png
