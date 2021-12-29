#!/bin/bash
linuxnd login -u $1 -p $2
mkdir data
linuxnd cp "oss://$3" -d data
