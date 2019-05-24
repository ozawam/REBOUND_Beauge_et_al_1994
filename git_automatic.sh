#!/bin/bash
MESSAGE=${1:-"add a subtle touch"}

git add problem.c analysis/ Makefile src/ start.sh git_automatic.sh  
git commit -m "${MESSAGE}"
#git push original master
git push origin hashes
