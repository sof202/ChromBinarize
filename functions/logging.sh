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

cat 1>&2 << ERROR_MESSAGE
[ERROR: ${current_time}]
${message}
ERROR_MESSAGE
}
