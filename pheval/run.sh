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
while getopts ":h:g:t:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      t) # Table Name
         TABLE=$OPTARG;;
      g) # HG version
         HG=$OPTARG;;
      p) # Phenotypic Version
         PHENOTYPEV=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
EXOMISER_VERSION=13.1.0
HG_VERSION=${HG:-hg19}
PHENOTYPE_VERSION=${PHENOTYPEV:-2209}
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
  --exomiser.data-directory=/exomiser-data/"$PHENOTYPE_VERSION"_phenotype/ \
  --exomiser.hg19.data-version=$PHENOTYPE_VERSION \
  --exomiser.phenotype.data-version=$PHENOTYPE_VERSION
}

## THIS PROCESS HAS BEEN DOING IN PYTHON
valid_table () {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    echo $LIST | tr "$DELIMITER" '\n' | grep -F -q -x -i "$VALUE"
}

TABLES="HP_HP_MAPPINGS HP_MP_MAPPINGS HP_ZP_MAPPINGS"


if [ -z "$TABLE" ]
then
  init_logs
  log "Running"
  setting_python_env
  pull_image
  exomiser_run
else
  if ! valid_table "$TABLES" " " $TABLE; then
    log "Invalid $TABLE"
    exit 1
  fi
  log "Second Run"
  log "Table" $TABLE
  exomiser_run
  log "Done"
fi
