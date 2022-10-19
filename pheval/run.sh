#!/usr/bin/env bash
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
EXOMISER_VERSION=13.1.0
HG_VERSION=hg19
DATA_VERSION=2209
DATA_FOLDER="$SCRIPT_DIR/../data/${DATA_VERSION}_phenotype"
EXOMISER_FOLDER="$SCRIPT_DIR/../lib/exomiser"

init_logs () {
  # Reset the log files
  printf '' > $SCRIPT_DIR/../info.log > $SCRIPT_DIR/../debug.log

  # Tail the info logfile as a background process so the contents of the
  # info logfile are output to stdout.
  tail -f $SCRIPT_DIR/../info.log &

  # Set an EXIT trap to ensure your background process is
  # cleaned-up when the script exits
  trap "pkill -P $$" EXIT
  # Redirect both stdout and stderr to write to the debug logfile
  exec 1>>debug.log 2> >(tee -a $SCRIPT_DIR/../debug.log >&2)
}

log () {
  # Write to both info and debug
  d=$(date "+%Y-%m-%d-%H:%M:%S")
  echo $d "- PHEVAL -" $1 | tee -a $SCRIPT_DIR/../info.log $SCRIPT_DIR/../debug.log
}

pull_image () {
  docker pull exomiser/exomiser-cli
}

download () {
  FILENAME=$(basename "$1")
  log "Downloading file $FILENAME"
  if test -f "$FILENAME"; then
    log "$FILENAME exists."
    return
  fi
  wget $1
  log "Unzipping $FILENAME"
  unzip $FILENAME
}

prepare_exomiser () {
  if test -d "$EXOMISER_FOLDER"; then
    log "$EXOMISER_FOLDER exists."
    log "Copying example files"
    cp -v $EXOMISER_FOLDER/examples/test-analysis-exome.yml $EXOMISER_FOLDER/exomiser-config/
    cp -v $EXOMISER_FOLDER/examples/Pfeiffer.vcf $EXOMISER_FOLDER/exomiser-config/
    return
  fi
  EXOMISER_FILE="exomiser-cli-$EXOMISER_VERSION-distribution.zip"
  download "https://data.monarchinitiative.org/exomiser/latest/$EXOMISER_FILE"
  unzip -u $EXOMISER_FILE -d ./
  mv ./exomiser-cli-$EXOMISER_VERSION/ $EXOMISER_FOLDER
  mkdir -p $EXOMISER_FOLDER/exomiser-config
  log "Copying example files"
  cp -v $EXOMISER_FOLDER/examples/test-analysis-exome.yml $EXOMISER_FOLDER/exomiser-config/
  cp -v $EXOMISER_FOLDER/examples/Pfeiffer.vcf $EXOMISER_FOLDER/exomiser-config/
}

prepare_data () {
  log "Preparing exomiser data"
  mkdir -p $DATA_FOLDER
  cd $DATA_FOLDER
  download "https://data.monarchinitiative.org/exomiser/latest/"$DATA_VERSION"_"$HG_VERSION".zip"
  download "https://data.monarchinitiative.org/exomiser/latest/"$DATA_VERSION"_phenotype.zip"
  cd -
}

setting_python_env () {
  log "Setting python environment"
  if test -d "$SCRIPT_DIR/../.venv"; then
    source $SCRIPT_DIR/../.venv/bin/activate
    return
  fi
  python -m venv venv
  source $SCRIPT_DIR/../.venv/bin/activate
  export PYTHONPATH=.:$PYTHONPATH
  pip install -e .
}

exomiser_run() {
  log "Running exomiser"
  docker run -v "$SCRIPT_DIR/../data/:/exomiser-data" \
  -v "$SCRIPT_DIR/../lib/exomiser/exomiser-config/:/exomiser"  \
  -v "$SCRIPT_DIR/../lib/exomiser/results:/results"  \
  exomiser/exomiser-cli:$EXOMISER_VERSION \
  --assembly $HG_VERSION \
  --analysis /exomiser/test-analysis-exome.yml  \
  --vcf /exomiser/Pfeiffer.vcf  \
  --exomiser.data-directory=/exomiser-data/"$DATA_VERSION"_phenotype/ \
  --exomiser.hg19.data-version=$DATA_VERSION \
  --exomiser.phenotype.data-version=$DATA_VERSION
}

init_logs
setting_python_env
prepare_exomiser
prepare_data
pull_image
exomiser_run
log "Done"