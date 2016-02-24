CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
van_leer

## Higher order extension: none or ppm
none

## CFL
0.5

## Time-stepping method: global (g) or local (l)
g

## Tolerance for residue norm comparisons
1e-6

## Grid file
bumpgrid-90-rot.txt

## State Load File ('~' for no load file)
###load.fvtk
~

## Max Iterations
25000

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
1000

## Debug level: Most detail (1) to least detail (5)
5

# FLOW SPECIFIC

## gamma (ratio of specific heats)
1.4

## R\_gas (specific gas constant)
287.

## FREE STREAM PROPERTIES

### Free Stream Density
0.00173

### Free Stream X Speed (at M = 2.4 and T = 300)
732.83

### Free Stream Y Speed
0.

### Free Stream Pressure
10.379


### Post Processing

### Plot type
Pseudocolor

### Plot variable
Mach_No
