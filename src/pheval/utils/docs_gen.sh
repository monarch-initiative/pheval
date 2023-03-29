#!/bin/bash

# set -e
cd ../../../
SOURCE_FOLDER='./src'
FILES=$(find $SOURCE_FOLDER -type f -iname '*.py' -not -iname '__init__.py' -not -empty)
rm -rf ./docs/api

for f in $FILES
do
    clean_dir=${f#./src/}
    last_folder=`dirname $clean_dir`
    full_fname="${f##*/}"
    fname="${full_fname%%.*}"
    mkdir -p ./docs/api/$last_folder
    ref=$(echo $f | sed 's#/#.#g' | sed 's/..src/src/g' | sed 's/\.[^.]*$//')
    echo ::: $ref >> ./docs/api/$last_folder/$fname.md
done