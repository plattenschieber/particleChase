#!/bin/bash

./particleChase > hans

gnuplot -persist << END_DATA
  set title 'Circle'
  set pointsize 0.00001
  plot 'hans' using 2:3 title 'particle 1' with p
END_DATA
