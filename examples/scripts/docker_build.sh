#!/bin/bash -e
# docker does not follow symlinks, so copy




# NOTE: basename will strip any username/ prefix, but this is intentional.
docker build -f "$(basename "$1")" -t "$(basename "$2")" .

# NOTE: The nested double quotes above are apparently correct.
# See https://unix.stackexchange.com/questions/100945/word-splitting-when-parameter-is-used-within-command-substitution