import numpy as np
import math

GAMMA = 1.14
gamma = GAMMA

G1 = (GAMMA - 1.0)/(2.0*GAMMA)
G2 = (GAMMA + 1.0)/(2.0*GAMMA)
G3 = 2.0*GAMMA/(GAMMA - 1.0)
G4 = 2.0/(GAMMA - 1.0)
G5 = 2.0/(GAMMA + 1.0)
G6 = (GAMMA - 1.0)/(GAMMA + 1.0)
G7 = (GAMMA - 1.0)/2.0
G8 = GAMMA - 1.0
G9 = (GAMMA + 1.0) / (GAMMA - 1.0)
G10 = GAMMA / (GAMMA - 1.0)
G11 = 1 / (GAMMA - 1.0)
G12 = (GAMMA - 1.0) / GAMMA

delta = (GAMMA - 1.0)/2.0

R_gas_const = 287.049

cp = GAMMA *R_gas_const / (GAMMA - 1)

