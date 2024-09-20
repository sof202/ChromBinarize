#!/bin/bash

logs() {
  # allows for some messages to be logged only if in debug mode
  is_logged=$1
  message=$2

  current_time=$(date +'%-d-%m ~ %T')

  if [[ is_logged -eq 1 ]]; then
cat << LOG_MESSAGE
[LOG: ${current_time}]
${message}

LOG_MESSAGE
  fi
}

errors() {
  message=$1
  current_time=$(date +'%-d-%m ~ %T')
  RED='[0;31m'
  NO_COLOUR='[0m'

cat 1>&2 << ERROR_MESSAGE
${RED}
[ERROR: ${current_time}]
${message}
${NO_COLOUR}
ERROR_MESSAGE
}
