import os, sys
from pyscipopt import Model

from contextlib import contextmanager

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

ZeroEpsilon = 10e-10

if __name__ == "__main__":
    tests = []
    with open('tests_miplib/tests.txt', 'r') as tests_f:
        for line in tests_f:
            if_tested, name, answer = line.split()
            if if_tested == "+":
                tests.append((name, float(answer)))
    for name, answer in tests:
        # solve original problem; write answer
        print(name)
        with suppress_stdout():
            milp = Model()
            milp.redirectOutput()
            problem_name = 'tests_miplib/' + name + '.mps'
            milp.readProblem(problem_name)
            milp.optimize()
            orig_answer = float(str(milp.getPrimalbound()))
            # generate cuts for problem
            os.system(f'python3 z-python-scripts/prepare_data.py {problem_name}')
            os.system('./Cutting_Plane_Generators files/data.txt')
            os.system(f'python3 z-python-scripts/write_cuts_in_lp.py files/cuts.txt {problem_name}')
            # sovle cutted problem; compare answer 
            cmilp = Model()
            cmilp.redirectOutput()
            cproblem_name = problem_name.replace(".mps", "_cutted.lp")
            cmilp.readProblem(cproblem_name)
            cmilp.optimize()
            cutted_answer = float(str(cmilp.getPrimalbound()))
        
        diff = abs(orig_answer - cutted_answer)
        if diff > ZeroEpsilon:
            print("Test " + name + " is not passed")    
            print(orig_answer, cutted_answer, diff)
        else:
            print("Test " + name + " is passed")
