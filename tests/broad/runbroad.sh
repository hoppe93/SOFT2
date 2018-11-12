#!/bin/bash
#
# Run the "broad" SOFT tests, conducting
# tests of how different SOFT modules
# work together.
# ######################################

##
# Run the test named '$name'.
#
# name: Name of test to run
# #######################
function runtest {
    case $1 in
        "orbits")
            orbits/run.sh
            ;;
        *)
            error "No test named '$1'."
            ;;
    esac
}

##
# Print a help message
# ###############
function hlp() {
    echo -e "Syntax: runbroad.sh test1 [test2 [test3 [...]]]\n"

    echo -e "If 'all' is provided as the sole argument, all tests are executed.\n"

    echo "AVAILABLE TESTS:"
    echo "  orbits          -- Test the orbits tool"
}

# MAIN PROGRAM
if [ "$#" -eq 1 ] && [ "$@" == "all" ]; then
    runall
elif [ "$#" -ge 1 ]; then
    for arg in "$@"
    do
        runtest $arg
    done
else
    hlp
    ./error "Hej"
    exit
fi

