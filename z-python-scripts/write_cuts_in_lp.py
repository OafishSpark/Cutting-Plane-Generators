from sys import argv
from os.path import dirname

from pyscipopt import Model

if __name__ == "__main__":
    try:
        _, input_path, task_path = argv
        print(input_path, task_path)
    except ValueError:
        input_path = "files/cuts.txt"
        task_path = "files/task.lp"
    model = Model()
    answer = ""
    model.readProblem(task_path)
    vars = model.getVars()
    with open(input_path, 'r') as file:
        i = 0
        for line in file:
            if line == "": break
            answer += "c" + str(i) + ": "
            i += 1
            elems = line.split()
            for elem in elems:
                if elem == ">=": 
                    answer = answer[:-2]
                    break
                ind, val = elem.split(",")
                answer += val + " " + str(vars[int(ind)]) + " + "
            answer += ">= " + elems[-1] + "\n"
    with open(task_path, 'r') as f:
        with open(task_path.replace(".lp", "_cutted.lp"), 'w') as out:
            for line in f:
                out.write(line)
                if line.strip() == "Subject To":
                    out.write(answer)
                
                
            

