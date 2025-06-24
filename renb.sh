#!/bin/bash
# This script opens a notebook file with nano editor

# Define the path to the notebook file
NOTEBOOK_FILE=~/Documents/notes/CPDprogress.docx 

# Open the notebook file with nano editor
wps "$NOTEBOOK_FILE" &
