#!/bin/bash

#this is the correct location in this workflow's container
export RATT_HOME=/venv/opt/RATT #remove when Dockerizeing

$RATT_HOME/start.ratt.sh $1 $2 $3 $4