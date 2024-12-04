# Quantum Error Correction 

This repository contains some introduction to quantum error correction (QEC), along with code for simulating various quantum error-correcting codes. Each subdirectory focuses on a specific code or code family with Jupyter notebooks and accompanying scripts.

## Project Structure

- **`3qubit_code/`**: Simulates the 3-qubit code for correcting simple bit flip errors, relaxation errors, and depolarization errors. Mainly for introductory purposes to QEC and quantum noise.
- **`Surface_Code/`**: Introduces the Surface Code using Qiskit, focusing on the [7-qubit error detection code](https://arxiv.org/abs/1912.09410).
- **`Bosonic_Modes/`**: A brief introduction to the square lattice [Gottesman-Kitaev-Preskill](https://errorcorrectionzoo.org/c/multimodegkp) (GKP) code, which encodes a qubit into an oscilaltor.
- **`repetition_cat_code/`**: Replicates numerical results for logical phase flip probabilities in the [repetition cat code](https://errorcorrectionzoo.org/c/cat_concatenated) (arxiv:2302.06639, Appendix E). The repetition cat code concatenates a cat qubit (encoded in an oscillator) with a repetition code.
- **`O2O/`**: Investigates [oscillator-to-oscillator](https://errorcorrectionzoo.org/c/oscillators_into_oscillators) (O2O) codes using Strawberry Fields. O2O codes aim to protect the entire Bosonic Hilbert space of the oscillator.

