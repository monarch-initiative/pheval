#!/usr/bin/env bash
set -e

Help()
{
   # Display Help
   echo "PhEval Runner"
   echo
   echo "Syntax: scriptTemplate [-t|h|g|p]"
   echo "options:"
   echo "t     Table name that will be scrambled  (valid options are: HP_HP_MAPPINGS; HP_MP_MAPPINGS; HP_ZP_MAPPINGS; ALL; - All option will scramble every mentioned table )"
   echo "g     Human Genome Version (valid options are: hg19; hg38 - Default version: hg19 )"
   echo "p     Human Genome Version (Default version: 2209 )"
   echo "h     Print this Help."
   echo
}

############################################################
############################################################
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":htg:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      t) # Table Name
         TABLE=$OPTARG;;
      g) # HG version
         HG=$OPTARG;;
      p) # Phenotypic Version
         DATAV=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
EXOMISER_VERSION=13.1.0
HG_VERSION=${HG:-hg19}
DATA_VERSION=${DATAV:-2209}
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

exomiser_run () {
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

dump () {
  log "Dumping $1"
  echo "CALL CSVWRITE('../output/$1.csv',  'SELECT * FROM EXOMISER.$1',  'charset=UTF-8 fieldSeparator=;');" > dump.sql
  java -Dh2.bindAddress=127.0.0.1 -cp "$(pwd)/../lib/h2.jar" org.h2.tools.RunScript -url jdbc:h2:file:$(pwd)/../data/2209_phenotype/2209_phenotype/2209_phenotype -script dump.sql -user sa
}

valid_table () {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    echo $LIST | tr "$DELIMITER" '\n' | grep -F -q -x -i "$VALUE"
}

table_dumping () {
  TABLES_ALL="all HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS"
  TABLES_ARRAY=( HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS )
  TABLES="HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS"

  if valid_table "$TABLES_ALL" " " $TABLE; then
      if [[ $(fgrep -ix $TABLE <<< "all")  ]]; then
        for i in "${TABLES_ARRAY[@]}"
        do
          dump $i  
        done
        return
      fi
      dump $TABLE
  else
      log "Invalid table - $TABLE"
      exit 1
  fi
}


init_logs
log "Running"
log "HG Version $HG_VERSION"
log "Phenotypic Data Version $DATA_VERSION"
if [ -z "$TABLE" ]
then
  log "First Run"
  setting_python_env
  prepare_exomiser
  prepare_data
  pull_image
  exomiser_run
else
  log "Second Run"
  log "Table" $TABLE
  table_dumping
  exomiser_run
  log "Done"
fi
