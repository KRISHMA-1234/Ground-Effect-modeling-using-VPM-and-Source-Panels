# Vortex Particle Code

## How to install

```bash
cd build
cmake ..
make
```

## Test the kernels

To test kernels with different order of expansion

```bash
./exafmm_kernel 10
./exafmm_kernel 20
./exafmm_kernel 30
```

## Run VPM simulation

```bash
cp ../exafmm_vpm/init.txt . # Copy init file to build directory
./exafmm # run exafmm
```

Output should be

```bash
--- FMM Profiling    ------------
Direct N-Body        : 0.052411 s
--- FMM vs. direct   ------------
Rel. L2 Error (Str)  : 5.65168e-01 s
Rel. L2 Error (Vel)  : 2.21746e-01 s
```
