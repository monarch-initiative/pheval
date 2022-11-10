# PhEval Contribution guide


## Contribution guidelines:

1. Every change must first be presented in an issue
1. Any change should be as small as possible (one "logical" unit, like adding a CLI method)
1. Every change is suggested in the form of a pull request and should be reviewed by at least 1 other PhEval developer. Changes to the CLI itself should be reviewed by at least 2 other PhEval developers.


## Repository structure

- All utility methods, i.e. methods that do not concern the actual PhEval runners (or other pipeline steps), should go in `src/pheval/utils.py`.

## Style guide

- There should be no functional code in the CLI method bodies (other then delegating to a library method)
- Methods involving IO (reading and writing files) should strictly seperate between the IO and processing part. I.e there should be a method that takes in the objects used for the processing _after_ reading them from wherever they are stored, returning the _object_ of whatever is being produced just before writing.
