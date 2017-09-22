#!/bin/bash

PIDFILE=/global/homes/e/erbilgin/repos/magi/magi_job/task.pid

if [ -f $PIDFILE ]
then
  PID=$(cat $PIDFILE)
  ps -p $PID > /dev/null 2>&1
  if [ $? -eq 0 ]
  then
    exit 1
  else
    ## Process not found assume not running
    echo $$ > $PIDFILE
    if [ $? -ne 0 ]
    then
      echo "Could not create PID file"
      exit 1
    fi
  fi
else
  echo $$ > $PIDFILE
  if [ $? -ne 0 ]
  then
    echo "Could not create PID file"
    exit 1
  fi
fi

/usr/common/software/python/2.7-anaconda-4.4/bin/python -W ignore /global/homes/e/erbilgin/repos/magi/magi_job/mirror_process_submit.py

rm $PIDFILE
