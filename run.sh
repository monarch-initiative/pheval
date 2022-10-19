#!/usr/bin/env bash

set -e

EXOMISER_VERSION=13.1.0
HG_VERSION=hg19
DATA_VERSION=2209
DATA_FOLDER="./data/${DATA_VERSION}_phenotype"
EXOMISER_FOLDER='./exomiser'

init_logs () {
    # Reset the log files
  printf '' > info.log > debug.log

  # Tail the info logfile as a background process so the contents of the
  # info logfile are output to stdout.
  tail -f info.log &

  # Set an EXIT trap to ensure your background process is
  # cleaned-up when the script exits
  trap "pkill -P $$" EXIT
  # Redirect both stdout and stderr to write to the debug logfile
  exec 1>>debug.log 2> >(tee -a debug.log >&2)
}

log () {
  # Write to both info and debug
  d=$(date "+%Y-%m-%d-%H:%M:%S")
  echo $d "- PHEVAL -" $1 | tee -a info.log debug.log
}

check_root () {
  if [[ $(whoami) != "root" ]]; then
    log 'Try to run it with sudo'
    exit 1
  fi
}

pull_image () {
  RUNNING_MSG=$( docker pull exomiser/exomiser-cli)
  SUB_RUNNING_MSG='Is the docker daemon running?'
  if [[ "$RUNNING_MSG" == *"$SUB_RUNNING_MSG"* ]]; then
    RUN_STATUS=$(service docker start)
    log "Starting docker"
    docker pull exomiser/exomiser-cli
  fi
  log 'Pulling exomiser docker image'
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
    cp -v ./exomiser/examples/test-analysis-exome.yml $EXOMISER_FOLDER/exomiser-config/
    cp -v ./exomiser/examples/Pfeiffer.vcf $EXOMISER_FOLDER/exomiser-config/
    return
  fi
  EXOMISER_FILE="exomiser-cli-$EXOMISER_VERSION-distribution.zip"
  download "https://data.monarchinitiative.org/exomiser/latest/$EXOMISER_FILE"
  unzip -u $EXOMISER_FILE -d ./
  mv ./exomiser-cli-$EXOMISER_VERSION/ $EXOMISER_FOLDER
  mkdir -p $EXOMISER_FOLDER/exomiser-config
  log "Copying example files"
  cp -v ./exomiser/examples/test-analysis-exome.yml $EXOMISER_FOLDER/exomiser-config/
  cp -v ./exomiser/examples/Pfeiffer.vcf $EXOMISER_FOLDER/exomiser-config/
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
  if test -d "./.venv"; then
    source ./.venv/bin/activate
    return
  fi
  python -m venv venv
  source $pwd/.venv/bin/activate
  export PYTHONPATH=.:$PYTHONPATH
  pip install -e .
}

exomiser_run() {
  log "Running exomiser"
  docker run -v "$(pwd)/data/:/exomiser-data" \
  -v "$(pwd)/exomiser/exomiser-config/:/exomiser"  \
  -v "$(pwd)/exomiser/results:/results"  \
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
check_root
pull_image
exomiser_run