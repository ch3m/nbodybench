import NbodyBench
import time
from colours import *
bench = NbodyBench.interface()
bench.setSimulation(10000,1000,0.003,100)
start_time = time.time()
bench.unitstride()
elapsed_time = time.time() - start_time
printout("elapsed_time: " + str(elapsed_time) + "\n", GREEN)
