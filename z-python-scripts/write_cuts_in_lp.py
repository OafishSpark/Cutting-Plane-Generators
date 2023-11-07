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
            answer += " c" + str(i) + ": "
            i += 1
            elems = line.split()
            for elem in elems:
                if elem == ">=": 
                    answer = answer[:-2]
                    break
                ind, val = elem.split(",")
                answer += val + " " + str(vars[int(ind)]) + " + "
            answer += ">= " + elems[-1] + "\n"
    task_path = task_path.replace(".lp", "_cutted.lp")
    task_path = task_path.replace(".mps", "_cutted.lp")
    with open(task_path, 'w') as fu:
        pass
    print(task_path)
    model.writeProblem(task_path)
    lines = []
    with open(task_path, 'r') as out:
        lines = out.readlines()
    with open(task_path, 'w') as out:
        for line in lines:
            out.write(line)
            if line == 'Subject to\n':
                out.write(answer)
