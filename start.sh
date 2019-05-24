#!/bin/bash

rm -r output/snap01
rm nohup.out

mkdir output
mkdir output/snap01

./rebound &
