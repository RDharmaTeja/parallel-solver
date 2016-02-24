CONFIG FILE
===========

# SOLVER SPECIFIC

## Scheme: van_leer / ldfss0 / hlle / ausm
van_leer

## Higher order extension: none or ppm or muscl
none

## CFL
0.5

## Time-stepping method: global (g) or local (l)
l

## Higher Order Time Integration: none or RK4
RK4

## Tolerance for residue norm comparisons
1e-6

## Grid file
bumpgrid-3D-1row-2.txt

## State Load File ('~' for no load file)
##state.fvtk
~


## Max Iterations
40000

## Checkpoint iter (dump data after how many iterations)
### (Enter 0 to turn checkpointing off)
1000

## Debug level: Most detail (1) to least detail (5)
1

# FLOW SPECIFIC

## gamma (ratio of specific heats)
1.4

## R\_gas (specific gas constant)
287.

## FREE STREAM PROPERTIES

### Number of variables
5

### Free Stream Density
0.1613

### Free Stream X Speed (at M = 2.4 and T = 300)
567.16

### Free Stream Y Speed
0.

### Free Stream Z Speed
0.

### Free Stream Pressure
5930.32

### Viscous effects
### Give mu reference = 0 for inviscid
### Using Sutherlands law for coefficient of viscosity
### mu reference or mu0 (in kg/ms)
0.02

### T reference or T0 (in K)
273.15

### Sutherland Temparature (K)
110.4

### Prandtl Number
0.7

### Post Processing

### Plot type
Pseudocolor

### Plot variable
Mach_No
