#!/bin/bash 
# Copy over file, run LIFT and tar result
scp -P 39578 $1 carl@129.16.26.50:~/Documents/Code/LIFT/LIFT/data/testimg/$1
ssh -p 39578 carl@129.16.26.50 "cd ~/Documents/Code/LIFT/LIFT; ./run.sh $1; tar -czvf results_$1.tar.gz results/"

# Copy result over and remove temporary files 
scp -P 39578 carl@129.16.26.50:~/Documents/Code/LIFT/LIFT/results_$1.tar.gz ./ 
ssh -p 39578 carl@129.16.26.50 "cd ~/Documents/Code/LIFT/LIFT; rm results_$1.tar.gz; cd results; rm *; cd ..; rm data/testimg/$1"


