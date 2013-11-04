import NbodyBench
import time
from colours import *
bench = NbodyBench.interface()
#bench.setSimulation(10000,5000,0.006,1000)
#configuracion[Numero de particulas, numero de pasos, deltaT, numero de promedios, temperatura]
bench.setSimulation(10,500,0.006,100,0.85)
bench.setBox(27.144178836455939,27.144178836455939,27.144178836455939)
bench.setBoxFilename("InitialBox")
bench.read_box()
bench.printParams()
start_time = time.time()
bench.unitstride()
elapsed_time = time.time() - start_time
printout("elapsed_time: " + str(elapsed_time) + "\n", GREEN)
if (elapsed_time <= 120 ):
  printout("Test (elapsed time) Passed...\n", RED)
else:
  printout("Test (elapsed time) Failed...\n", BLUE)
