import numpy as np

mu = 0.082

V = 67

d = 0.075

rho = 1.2

Re = rho * V * d / mu

print 'Re = ', Re
print 'T = ', 98538/(rho * 286.6)

x = 0.04

Re_x = rho * V * x / mu
blt = 4.91 * x / np.sqrt(Re_x)

print 'BLT = ', blt
