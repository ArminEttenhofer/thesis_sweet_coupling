#! /bin/bash

#
# This is a helper script for development environments.
#
# It loads the SWEET environment and then starts scons with the parameters provided.
#

# Location of this script
SCRIPTDIR="$(dirname "${BASH_SOURCE[0]}")"
cd "${SCRIPTDIR}"
cd ".."


# We need to put this into an env to ensure that the program arguments are not forwarded
activateenv() {
	source "activate.sh"
}

activateenv

scons $@
