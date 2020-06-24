#!/usr/bin/env bash

#This is a helper script to install CCP4, configure the binaries and add it to the path.
#CCP4, while open sourced and mostly GNU licensed (certainly, the tool STACEI uses are at least)
#please take some time to refer to the CCP4 license at: https://www.ccp4.ac.uk/commercial/main.php


url="ftp://ftp.ccp4.ac.uk/ccp4/7.0/ccp4-7.0-linux-x86_64.tar.bz2"
file="ccp4-7.0-linux-x86_64.tar.bz2"
dir="ccp4-7.0"

#download, unzip and enter ccp4
echo "Downloading CCP4"
wget "$url"

echo "Untarring it"
tar -xjf "$file"

echo "Copying CCP4 -> stacei_src/bin"
if [ -d "stacei_src/bin/$dir" ]; then
  echo "Deleteing old directory"
  rm -rf "stacei_src/bin/$dir"
fi

echo "Moving to STACEI bin"
mv "$dir" stacei_src/bin

#run binary setup
echo "Running binary setup"
bash stacei_src/bin/"$dir"/BINARY.setup

echo "Removing tarbell"
rm "$file" -rf