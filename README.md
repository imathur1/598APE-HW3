# 598APE-HW3

This repository contains code for homework 3 of 598APE.

This assignment is relatively simple in comparison to HW1 and HW2 to ensure you have enough time to work on the course project.

In particular, this repository is an implementation of an n-body simulator.

To compile the program run:
```bash
make -j
```

To clean existing build artifacts run:
```bash
make clean
```

To compile the program with profiling on run:
```bash
make -j PROFILE=1
```

To run the tests to check for correctness:
```bash
make test
./run_tests.sh
```

If compiled with profiling, running the program will output a `my_profile.prof` file which can be visualized by running:
```bash
make view-profile
```


This program assumes the following are installed on your machine:
* A working C compiler (g++ is assumed in the Makefile)
* make

The nbody program is a classic physics simulation whose exact results are unable to be solved for exactly through integration.

Here we implement a simple time evolution where each iteration advances the simulation one unit of time, according to Newton's law of gravitation.

Once compiled, one can call the nbody program as follows, where nplanets is the number of randomly generated planets for the simulation, and timesteps denotes how long to run the simulation for:
```bash
./main.exe <nplanets> <timesteps>
```

In particular, consider speeding up simple run like the following (which runs ~6 seconds on my local laptop under the default setup):
```bash
./main.exe 1000 5000
```

Exact bitwise reproducibility is not required, but approximate correctness (within a reasonable region of the final location).

# Optimizations
A full analysis can be found in our paper, `Mini_Paper_3__Nbody.pdf`. Our charts were generated in analysis.py using the following commands:
1. Install pip:
   ```bash
   sudo apt install python3-pip
   ```
2. Install matplotlib:
    ```bash
    pip install matplotlib
    ```
3. Run the analysis:
    ```bash
    python analysis.py
    ```

You can run our optimizations with the following commands:
1. Checkout commit hash:
   ```bash
   git checkout <hash>
   ```
2. Recompile program:
    ```bash
    make clean && make
    ```
3. Run simulation:
    ```bash
    ./main.exe 1000 5000
    ```
