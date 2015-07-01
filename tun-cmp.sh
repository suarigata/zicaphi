#!/bin/bash
# /***************************************************************************
#  *
#  * Copyright 2015 Ian Liu Rodrigues and Edson Borin
#  *
#  ***************************************************************************/
function fail
{
    echo "ERROR: $@"
    exit 1
}

TIME=${TIME:-`which time`}
DATA=${DATA:-/home/edson/cmp/sn3d.su}
ARGS=${ARGS:-"1.98e-7 1.77e-6 101 600 0.002"}
PROG=${PROG:-./cmp.x}

# Check time program
[ -f "$PROG" ] || fail "$PROG is not a valid program (Did you build it?)."

# Check cmp program
[ -f "$TIME" ] || fail "$TIME is not a valid program."

# Check input
[ -f "$DATA" ] || fail "$DATA is not a valid input file (Did you copy it here?)."

# Run it
$TIME $PROG $ARGS $DATA
