#!/bin/bash
# should be called in the project directory
rm ./b-inv.txt
python3 ./b-inv-part/main.py ./task.lp >> ./b-inv.txt
rm ./task.mps
#./cuts.exe
