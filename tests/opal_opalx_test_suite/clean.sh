#!/bin/bash

# Set the path to the folder you want to clean

echo "Cleaning folder: 'opalx'"
echo "Deleting files that *do not* have a .in extension..."

# Use find to locate the files and delete them
rm -v opalx/*
rm -v *.png

echo "Cleaning complete."