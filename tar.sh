#!/bin/bash
 
if [[ $1 = '' ]]; then
    echo "missing input 'tarpv inputfile/inputfold'"
    exit 1
fi
 
tar -cf - "$1" | pv -s $(($(du -sk "$1" | awk '{print $1}') * 1024)) | gzip > $1.tar.gz
