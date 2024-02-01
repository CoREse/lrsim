#!/usr/bin/bash
for modelfile in models/*.gz; do
    gunzip -c $modelfile > models/`basename $modelfile .gz`
done