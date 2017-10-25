#!/bin/sh

if [ "$1" = "-show" ]; then
  ctest -N
else
  polecenie="ctest -V -S ../test.cmake"
  for nazwa in $*
  do
    polecenie="$polecenie,$nazwa"
  done
  echo $polecenie
  $polecenie
fi

