#!/bin/bash
# should be called in the project directory
rm ./files/b-inv.txt
python3 ./b-inv-part/main.py ./files/task.lp >> ./files/b-inv.txt
rm ./files/task.mps
#./cuts.exe
