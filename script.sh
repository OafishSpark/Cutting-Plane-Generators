#!/bin/bash
# should be called in the project directory
conda activate cutting-planes
python3 ./z-python-scripts/prepare_data.py "./files/plane.lp"
#./build/Cutting_Plane_Generators files/data.txt
#python3 ./z-python-scripts/write_cuts_in_lp.py ./files/cuts.txt ./files/cutted_task.lp
