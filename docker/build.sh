#!/bin/bash
cp -R ../workflow . &&
rm -R workflow/*.zip &&
cp -R ../local_settings . &&
cp -R ../tests . &&
cp -R ../__init__.py . &&
cp blastp workflow/blastbin &&
cp makeblastdb workflow/blastbin &&
cp ../magi_env.yml .

echo building...
docker build -t magi .

# reset
rm -r workflow
rm -r local_settings
rm -r tests
rm -r __init__.py
rm magi_env.yml