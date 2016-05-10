#!/bin/bash
gcc -Wall -Wextra -pedantic -std=c11  argmc.c -o execargmc.out -lm
#valgrind -v ./execargmc.out
#./execargmc.out
