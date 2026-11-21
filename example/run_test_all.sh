#!/bin/bash
RED='\033[0;31m'
BLUE='\033[0;34m'
GREEN='\033[0;32m'
NC='\033[0m'

printf "${RED}This script runs for a test calculation. If you want to get a reliable result, more data are required.${NC} \n"
printf "${RED} You should activate conda environment first. If you use a dedicated Julia project environment, set 
  julia_interpreter = julia --project=JULIA/PROJECT/PATH
in the example configuration files. Replace JULIA/PROJECT/_PATH} with the actual path to your Julia project environment.${NC} \n"

printf "${BLUE}Example data are stored at dft ${NC} \n"
printf "${BLUE} tar -zxvf dft.tar.gz ${NC} \n"
tar -zxvf dft.tar.gz

printf "${BLUE} cd preprocess ${NC} \n"
cd "preprocess"

printf "${GREEN} 1. deeph-preprocess --config  preprocess.ini  ${NC} \n"
deeph-preprocess --config  preprocess.ini
printf "${BLUE} cd ../train ${NC} \n"
cd ../train

printf "${GREEN} 2. deeph-train --config  train.ini  ${NC} \n"
deeph-train --config  train.ini 

printf "${BLUE} cd ../inference ${NC} \n"
cd ../inference 

printf "${BLUE} tar -zxvf OLP.tar.gz ${NC} \n"
tar -zxvf OLP.tar.gz

printf "${GREEN} 3. inference-jx --config  inference.ini  ${NC} \n"
inference-jx --config  inference.ini
