RESOURCE='resources/Makefile.j2'
CONFIG='resources/pheval-config.yaml'
OUTPUT='Makefile'

Help()
{
   # Display Help
   echo "Utility for Makefile generation"
   echo
   echo "Syntax: generatemakefile.sh [--resource|--config|--output]"
   echo "options:"
   echo "--resource Makefile.j2 base file - Default: resources/Makefile.j2"
   echo "--config config.yaml file that will be used to file Makefile.j2. - Default: resources/pheval-config.yaml"
   echo "--output where the new makefile will be written - Default: Makefile"
   echo
}

while [ $# -gt 0 ]; do
  case "$1" in
    --resource=*)
      RESOURCE="${1#*=}"
      ;;
    --config=*)
      CONFIG="${1#*=}"
      ;;
    --output=*)
      OUTPUT="${1#*=}"
      ;;
    *)
      Help
      exit 1
  esac
  shift
done


j2 $RESOURCE $CONFIG > $OUTPUT

echo "Done"