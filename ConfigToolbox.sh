#!/bin/bash 
# ConfigToolbox.sh
# Run this script ONCE (when you first clone the git repository) to configure local settings

# =======================================================================  CONFIG DIRECTORIES
read -p "Enter directory to save results:" resultsPath
echo $resultsPath > code/SimulationResults.path

read -p "Enter directory of installed Lightspeed toolbox:" lsLibPath
echo $lsLibPath > code/LightspeedLibrary.path



