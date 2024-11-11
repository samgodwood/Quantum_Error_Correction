# Useful functions

import numpy as np
from scipy.linalg import sqrtm

def tensor_product(*args):
    # Assume args is a list of matrices
    result = np.eye(1)
    for mat in args:
        result = np.kron(result, mat)
    return result

# Define basis states
ket_0 = np.array([[1], [0]])
ket_1 = np.array([[0], [1]])

# Define Pauli gates
X_gate = np.array([[0, 1], [1, 0]])
Y_gate = np.array([[0,1j],[-1j,0]])
Z_gate = np.array([[1,0],[0,-1]])
Identity = np.identity(2)

#Extra gates
Hadamard = 1/np.sqrt(2)*np.array([[1,1],[1,-1]])
                               
#Define CNOT gate
cnot_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])

# Tensor product to create a 3-qubit CNOT with control on qubit 1 and target on qubit 2
# Tensor product to create a 3-qubit CNOT with control on qubit 2 and target on qubit 3
CNOT_12 = tensor_product(cnot_matrix, Identity)
CNOT_23 = tensor_product(Identity, cnot_matrix)

#5 qubit CNOT definitions

#First, qubit 0 (control) qubit 3 (targer)
# Create a 16x16 numpy array of zeros for the CNOT operation
cnot_array = np.zeros((16, 16), dtype=int)

# Mapping from CNOT operation
cnot_mappings = {
    16: 18, 17: 19, 18: 16, 19: 17,
    20: 22, 21: 23, 22: 20, 23: 21,
    24: 26, 25: 27, 26: 24, 27: 25,
    28: 30, 29: 31, 30: 28, 31: 29
}

# Set the 'target' elements to 1 in the CNOT array
for source, target in cnot_mappings.items():
    cnot_array[source-16, target-16] = 1  # Adjusting for zero-based indexing

# Create a 32x32 identity matrix
CNOT_0_3    = np.eye(32, dtype=int)

# Replace the lower 16x16 block of i32 with the CNOT array
CNOT_0_3[16:32, 16:32] = cnot_array

#--------------------------------------------------------------------------------------------------------------------------------------------

# qubit 1 (control) qubit 3 (target)

# Mapping for CNOT operation with qubit 1 as control and qubit 3 as target
cnot_1_3_mappings = {
    8: 10, 9: 11, 10: 8, 11: 9,
    12: 14, 13: 15, 14: 12, 15: 13,
    24: 26, 25: 27, 26: 24, 27: 25,
    28: 30, 29: 31, 30: 28, 31: 29
}

# Create a 32x32 identity matrix
cnot_1_3_identity_array = np.eye(32, dtype=int)

# Update only the specific elements based on the CNOT mapping
for source, target in cnot_1_3_mappings.items():
    # Set the original identity element to 0
    cnot_1_3_identity_array[source, source] = 0
    # Set the new target element to 1
    cnot_1_3_identity_array[source, target] = 1

CNOT_1_3 = cnot_1_3_identity_array

#--------------------------------------------------------------------------------------------------------------------------------------------

#qubit 1 (control) qubit 4 (target)

# Mapping for CNOT operation with qubit 1 as control and qubit 4 as target (0 through 31)
cnot_1_4_mappings = {
    8: 9, 9: 8, 10: 11, 11: 10,
    12: 13, 13: 12, 14: 15, 15: 14,
    24: 25, 25: 24, 26: 27, 27: 26,
    28: 29, 29: 28, 30: 31, 31: 30
}

# Create a 32x32 identity matrix
cnot_1_4_identity_array = np.eye(32, dtype=int)

# Update only the specific elements based on the CNOT mapping
for source, target in cnot_1_4_mappings.items():
    # Set the original identity element to 0 for affected states
    cnot_1_4_identity_array[source, source] = 0
    # Set the new target element to 1
    cnot_1_4_identity_array[source, target] = 1

CNOT_1_4 = cnot_1_4_identity_array


#--------------------------------------------------------------------------------------------------------------------------------------------

#qubit 2 (control) qubit 4 (target)

# Mapping for CNOT operation with qubit 1 as control and qubit 4 as target (0 through 31)
cnot_2_4_mappings = {
    4:5, 5:4, 6:7, 7:6,
    12:13, 13:12, 14:15, 15:14,
    20:21, 21:20, 22:23, 23:22,
    28:29, 29:28, 30:31, 31:30
    
}

# Create a 32x32 identity matrix
cnot_2_4_identity_array = np.eye(32, dtype=int)

# Update only the specific elements based on the CNOT mapping
for source, target in cnot_2_4_mappings.items():
    # Set the original identity element to 0 for affected states
    cnot_2_4_identity_array[source, source] = 0
    # Set the new target element to 1
    cnot_2_4_identity_array[source, target] = 1

CNOT_2_4 = cnot_2_4_identity_array

def apply_error(state, operation):
    """Applies a given error operation to the state."""
    return np.dot(operation, state)

def prep_error_channel(state, kappa_ratio, alpha_squared):
    """
    State preparation error channel where each qubit independently has a probability p of a Z error.
    
    Args:
        state (array): The state.
        kappa_ratio (float): k1/k2 (k1 is the single excitation loss rate, k2 is the engineered two-excitation dissipation rate).
        alpha_squared (float): Average cat qubit photon number.

    Returns:
        array: The modified state after potential error application.
    """
    # Probability of Z error on an individual qubit
    p = alpha_squared * kappa_ratio

    # Define the I and Z gates
    I = Identity  # Assuming Identity is predefined in your code
    Z = Z_gate    # Assuming Z_gate is predefined as well
    
    # List of possible operations
    operations = [I, Z]
    
    # Determine errors independently for each qubit
    error_indices = [np.random.choice([0, 1], p=[1-p, p]) for _ in range(3)]
    
    # Select the error operators based on random choice
    error_operators = [operations[i] for i in error_indices]
    
    # Construct the full error operator for the 3-qubit state
    error_to_apply = tensor_product(*error_operators)
    
    # Apply the error operator to the state
    error_state = np.dot(error_to_apply, state)
    
    return error_state


def cnot_error_channel(state, kappa_ratio, alpha_squared, control_index, target_index):
    """
    CNOT error channel applied after each CNOT gate in the circuit.
    Introduces single-qubit Z errors on the control and target qubits
    with probability given in Table I of arxiv:2302.06639

    Args:
        state (array): The full quantum state vector.
        kappa_ratio (float): k1/k2 (k1 is the single excitation loss rate, k2 is the engineered two-excitation dissipation rate).
        alpha_squared (float): Average cat qubit photon number.
        control_index (int): Index of the control qubit in the state vector.
        target_index (int): Index of the target qubit in the state vector.

    Returns:
        array: The modified state after potential CNOT-induced errors.
    """
    # Calculate probabilities for Z errors
    p_Z1 = alpha_squared * kappa_ratio + (np.pi**2) / (64 * alpha_squared)
    p_Z2 = 0.5 * alpha_squared * kappa_ratio

    # Define the identity and Z gates
    I = np.array([[1, 0],
                  [0, 1]], dtype=complex)
    Z = np.array([[1, 0],
                  [0, -1]], dtype=complex)

    # Randomly decide if a Z error occurs on the control qubit
    error_index_control = np.random.choice([0, 1], p=[1 - p_Z1, p_Z1])
    # Randomly decide if a Z error occurs on the target qubit
    error_index_target = np.random.choice([0, 1], p=[1 - p_Z2, p_Z2])

    # Possible operations for the error channels on each qubit
    operations = [I, Z]

    # Initialize the list of operators for each qubit
    num_qubits = int(np.log2(len(state)))
    error_operators = [I] * num_qubits  # Start with Identity for all qubits

    # Apply the error operator to the control qubit
    error_operators[control_index] = operations[error_index_control]
    # Apply the error operator to the target qubit
    error_operators[target_index] = operations[error_index_target]

    # Construct the full error operator as tensor product
    error_to_apply = error_operators[0]
    for op in error_operators[1:]:
        error_to_apply = np.kron(error_to_apply, op)

    # Apply the error operator to the state
    error_state = np.dot(error_to_apply, state)

    return error_state

def encode_logical_state(input_qubit, kappa_ratio, alpha_squared):
    """
    Encodes the input qubit into a logical state using the phase flip code.

    Args:
        input_qubit (array): State to encode with the repetition code.
        kappa_ratio (float): k1/k2 (k1 is the single excitation loss rate, k2 is the engineered two-excitation dissipation rate).
        alpha_squared (float): Average cat qubit photon number.

    Returns:
        array: The logical encoded quantum state.
    """
    # Initial three-qubit state with the input qubit and two ancillary |0⟩ qubits
    three_qubit_state = tensor_product(input_qubit, ket_0, ket_0)

    # Apply state preparation error channel
    three_qubit_state = prep_error_channel(three_qubit_state, kappa_ratio, alpha_squared)

    # Apply the first CNOT gate between qubits 0 and 1
    entangled_state = np.dot(CNOT_12, three_qubit_state)

    # Apply CNOT error channel after the first CNOT gate (qubits 0-1)
    entangled_state = cnot_error_channel(entangled_state, kappa_ratio, alpha_squared, control_index=0, target_index=1)

    # Apply the second CNOT gate between qubits 1 and 2
    entangled_state = np.dot(CNOT_23, entangled_state)

    # Apply CNOT error channel after the second CNOT gate (qubits 1-2)
    entangled_state = cnot_error_channel(entangled_state, kappa_ratio, alpha_squared, control_index=1, target_index=2)

    # Apply Hadamard gates to complete the logical encoding
    hadamard_all = tensor_product(Hadamard, Hadamard, Hadamard)
    logical_encoded_state = np.dot(hadamard_all, entangled_state)

    return logical_encoded_state


def add_ancilla_qubits(error_state, kappa_ratio, alpha_squared):
    """
    Adds two ancilla qubits to be used for stabilizer measurements

    Parameters:
        error_state (np.array): The current state of the logical qubits to which ancilla qubits will be entangled with.
        kappa_ratio (float): The ratio k1/k2, where k1 is the single excitation loss rate, and k2 is the engineered 
                             two-excitation dissipation rate.
        alpha_squared (float): The average photon number for the cat qubit, a parameter controlling error rates.

    Returns:
        np.array: The final state.
    """
    # ading two ancilla qubits initialized to |0⟩
    extended_state = tensor_product(error_state, ket_0, ket_0)
    
    # Apply Hadamard gates to change basis
    hadamard_on_data = tensor_product(Hadamard, Hadamard, Hadamard, Identity, Identity)
    pre_entangled_state = np.dot(hadamard_on_data, extended_state)
    
    # Apply CNOT gates between data and ancilla qubits for stabilizer measurements
    # CNOT_0_3: control on qubit 0, target on ancilla qubit 3
    # CNOT_1_3: control on qubit 1, target on ancilla qubit 3
    # CNOT_1_4: control on qubit 1, target on ancilla qubit 4
    # CNOT_2_4: control on qubit 2, target on ancilla qubit 4

    entangled_state = np.dot(CNOT_0_3, pre_entangled_state)
    entangled_state = cnot_error_channel(entangled_state, kappa_ratio, alpha_squared, control_index=0, target_index=3)

    entangled_state = np.dot(CNOT_1_3, entangled_state)
    entangled_state = cnot_error_channel(entangled_state, kappa_ratio, alpha_squared, control_index=1, target_index=3)

    entangled_state = np.dot(CNOT_1_4, entangled_state)
    entangled_state = cnot_error_channel(entangled_state, kappa_ratio, alpha_squared, control_index=1, target_index=4)

    entangled_state = np.dot(CNOT_2_4, entangled_state)
    entangled_state = cnot_error_channel(entangled_state, kappa_ratio, alpha_squared, control_index=2, target_index=4)

    # Revert back to X basis
    final_state = np.dot(hadamard_on_data, entangled_state)
    
    return final_state

def measure_rightmost_2_qubits(state):
    """
    Measures the rightmost two qubits of a 5 qubit state.

    Parameters:
        state (np.array): State to measure.

    Returns:
        str: The measurement outcome for the rightmost two qubits as a string ('00', '01', '10', or '11').
    """  

    # Define projection operators for two qubits
    P_00 = tensor_product(ket_0, ket_0).dot(tensor_product(ket_0, ket_0).T)
    P_01 = tensor_product(ket_0, ket_1).dot(tensor_product(ket_0, ket_1).T)
    P_10 = tensor_product(ket_1, ket_0).dot(tensor_product(ket_1, ket_0).T)
    P_11 = tensor_product(ket_1, ket_1).dot(tensor_product(ket_1, ket_1).T)
    
    # Tensor projection operators with identity for the 3 leftmost qubits
    I_3 = np.eye(8)  # Identity for 3 qubits
    P_00_full = np.kron(I_3, P_00)
    P_01_full = np.kron(I_3, P_01)
    P_10_full = np.kron(I_3, P_10)
    P_11_full = np.kron(I_3, P_11)
    
    # Compute probabilities for each outcome
    p_00 = np.abs(state.T @ P_00_full @ state)[0, 0] #[0,0] to extract just a number
    p_01 = np.abs(state.T @ P_01_full @ state)[0, 0]
    p_10 = np.abs(state.T @ P_10_full @ state)[0, 0]
    p_11 = np.abs(state.T @ P_11_full @ state)[0, 0]

    # Ensure probabilities sum up to 1 (handling minor numerical inaccuracies)
    total_prob = p_00 + p_01 + p_10 + p_11
    p_00 /= total_prob
    p_01 /= total_prob
    p_10 /= total_prob
    p_11 /= total_prob

    # Randomly choose a measurement outcome based on the probabilities
    outcomes = ['00', '01', '10', '11']
    result = np.random.choice(outcomes, p=[p_00, p_01, p_10, p_11])

    return result


import pymatching
def mwpm_correction(results_list, state):
    """
    Corrects a given quantum state based on syndrome results from a 3-qubit phase-flip code using minimum-weight
    perfect matching (MWPM).
    
    Args:
        results_list (list): List of binary strings representing syndromes from measurement cycles.
        state (np.array): Quantum state vector to apply corrections to.
        
    Returns:
        np.array: Corrected quantum state.
    """
    # Define check matrix for 3-qubit phase-flip code and initialize the matcher
    H_z = np.array([[1, 1, 0], [0, 1, 1]])
    zmatching = pymatching.Matching(H_z)
    
    # Initialize cumulative error tracker
    cumulative_errors = np.zeros(3, dtype=int)
    # Decode each cycle’s syndrome and update cumulative error pattern
    for binary_str in results_list:
        Z_syndrome = [int(bit) for bit in binary_str]
        zprediction = zmatching.decode(Z_syndrome)
        cumulative_errors = (cumulative_errors + zprediction) % 2

    # Define Pauli-Z corrections for each qubit
    Z_gate = np.array([[1, 0], [0, -1]])  # Pauli-Z gate
    ZII = np.kron(Z_gate, np.eye(4))       # Z on qubit 1
    IZI = np.kron(np.eye(2), np.kron(Z_gate, np.eye(2)))  # Z on qubit 2
    IIZ = np.kron(np.eye(4), Z_gate)       # Z on qubit 3
    
    # Apply corrections based on cumulative error pattern
    if cumulative_errors[0] == 1:
        state = np.dot(ZII, state)
    if cumulative_errors[1] == 1:
        state = np.dot(IZI, state)
    if cumulative_errors[2] == 1:
        state = np.dot(IIZ, state)

    return state
    

def discard_ancilla_from_statevector(state_vector):
    """Discard the ancilla qubits from the statevector and return the state of the 3 logical qubits."""
    tensor_representation = state_vector.reshape([2]*5)  # Reshape to a 5-qubit tensor
    reduced_state = np.sum(tensor_representation, axis=(3,4))  # Sum over the ancilla qubits
    flattened_state = reduced_state.flatten().reshape((8,1))  # Flatten and reshape to column vector
    return flattened_state/np.linalg.norm(flattened_state)  # Normalize and return
