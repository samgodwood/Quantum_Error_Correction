# Quantum Error Correction Simulation

This repository contains code and resources for simulating quantum error-correcting codes.

## Project Structure

- `Basic 3 qubit code.ipynb`: Jupyter notebook containing the simulation code for the 3-qubit bit flip code for simple bit flip errors. The simulation is carried out using linear algebra with numpy, and the results are visualized to understand the performance of the error correction code against random bit flip errors.

- `3 Qubit Code with Relaxation Errors.ipynb`: Jupyter notebook building upon the functions defined in `Basic 3 qubit code.ipynb`. This notebook focuses on simulating relaxation errors (T_1 errors) and plots the lifetime of physical qubits vs. logical qubits. It provides insights into how the error correction code performs under relaxation error scenarios.
