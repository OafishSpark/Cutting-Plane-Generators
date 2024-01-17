import gurobipy as gp
from gurobipy import GRB

m = gp.read("/mnt/d/dev/experiment101/model.mps")

# Solve

# m.setParam("Cuts", 0)
# m.setParam("MIRCuts", 0)
# m.setParam("RelaxLiftCuts", 0)
# m.setParam("RLTCuts", 0)
# m.setParam("GomoryPasses", 0)
# m.setParam("FlowCoverCuts", 0)
# m.setParam("ImpliedCuts", 0)

# 1. Flow-Cover (94)
# 2. Gomory (102/47/+RLT=root)
# 3. Implied Bounds (131/61/45)
# 4. MIR (126/93/55/root)

m.optimize()
m.printStats()
