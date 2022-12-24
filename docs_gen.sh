#!/bin/bash

# set -e

SOURCE_FOLDER='./src'
FILES=$(find $SOURCE_FOLDER -type f -iname '*.py' -not -iname '__init__.py' -not -empty)
rm -rfv ./docs/api

for f in $FILES
do

    clean_dir=${f#./src/}
    folder=${clean_dir%/*}
    ref=$(echo $f | sed 's#/#.#g' | sed 's/..src/src/g' | sed 's/\.[^.]*$//')
    mkdir -p ./docs/api/$folder/
    touch ./docs/api/$folder/reference.md
    echo  ::: $ref >> ./docs/api/$folder/reference.md
done