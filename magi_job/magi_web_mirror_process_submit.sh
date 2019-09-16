#!/bin/bash

PIDFILE='/global/homes/p/pasteur/repos/magi/magi_job/task.pid'

# does PIDFILE exist as a file?
if [ -f $PIDFILE ]
then
  PID=$(cat $PIDFILE)
  # check if pid is actually running
  ps -p $PID > /dev/null 2>&1
  # if it is running, exit
  if [ $? -eq 0 ]
  then
    exit 1
  else
    ## PIDFILE exists but not running
    # write PID to PIDFILE
    echo $$ > $PIDFILE
    # if echo fails, tell me about it
    if [ $? -ne 0 ]
    then
      echo "Could not create PID file"
      exit 1
    fi
  fi
else
  # write PID to PIDFILE
  echo $$ > $PIDFILE
  # if echo fails, tell me about it
  if [ $? -ne 0 ]
  then
    echo "Could not create PID file"
    exit 1
  fi
fi

/global/common/software/m2650/python-cori/bin/python -W ignore /global/u2/p/pasteur/repos/magi/magi_job/mirror_process_submit.py
#/usr/common/software/python/2.7-anaconda-4.4/bin/python -W ignore /global/u2/p/pasteur/repos/magi/magi_job/mirror_process_submit.py
rm $PIDFILE
