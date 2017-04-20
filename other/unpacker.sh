#!/bin/bash

i=0

while [ -f "$i.tar.gz" ]; do
    tar -xzf $i.tar.gz
    let i=i+1
done
