import numpy as np
from scipy.linalg import sqrtm

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

#--------------------------------------------------------------------------------------------------------------------------------------------

def tensor_product(*args):
    """
    Compute the tensor product of an arbitrary number of matrices.

    Parameters:
    *args: An unpacked sequence of np.ndarray matrices.

    Returns:
    np.ndarray: The resulting matrix after computing the tensor product.
    """
    # Initialize with the identity matrix of size 1x1
    result = np.eye(1)
    # Compute the Kronecker product for each matrix in args
    for mat in args:
        result = np.kron(result, mat)
    return result

#--------------------------------------------------------------------------------------------------------------------------------------------

#Define CNOT gate
cnot_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])

#--------------------------------------------------------------------------------------------------------------------------------------------

#5 qubit CNOT definitions

# Create a 3-qubit CNOT gate with control on qubit 1 and target on qubit 2
CNOT_0_1 = tensor_product(cnot_matrix, Identity)

# Create a 3-qubit CNOT gate with control on qubit 2 and target on qubit 3
CNOT_1_2 = tensor_product(Identity, cnot_matrix)

#--------------------------------------------------------------------------------------------------------------------------------------------
 
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

#--------------------------------------------------------------------------------------------------------------------------------------------

#Define depolarizing error
def depolarizing_error(rho, p, qubit=None, all=False):
    k_0 = np.sqrt(1 - p) * Identity
    k_1 = np.sqrt(p/3) * X_gate
    k_2 = np.sqrt(p/3) * Z_gate
    k_3 = np.sqrt(p/3) * Y_gate

    num_qubits = int(np.log2(rho.shape[0]))

    if all:
        for i in range(num_qubits):
            K0 = tensor_product(np.eye(2**i), k_0, np.eye(2**(num_qubits - i - 1)))
            K1 = tensor_product(np.eye(2**i), k_1, np.eye(2**(num_qubits - i - 1)))
            K2 = tensor_product(np.eye(2**i), k_2, np.eye(2**(num_qubits - i - 1)))
            K3 = tensor_product(np.eye(2**i), k_3, np.eye(2**(num_qubits - i - 1)))
        
            rho = K0 @ rho @ K0.T + K1 @ rho @ K1.T + K2 @ rho @ K2.T + K3 @ rho @ K3.conj().T
    
        return rho

    else:
        if qubit is None:
            raise ValueError("qubit must be specified when applying to specific qubit")

        K0 = tensor_product(np.eye(2**qubit), k_0, np.eye(2**(num_qubits - qubit - 1)))
        K1 = tensor_product(np.eye(2**qubit), k_1, np.eye(2**(num_qubits - qubit - 1)))
        K2 = tensor_product(np.eye(2**qubit), k_2, np.eye(2**(num_qubits - qubit - 1)))
        K3 = tensor_product(np.eye(2**qubit), k_3, np.eye(2**(num_qubits - qubit - 1)))

        rho = K0 @ rho @ K0.T + K1 @ rho @ K1.T + K2 @ rho @ K2.T + K3 @ rho @ K3.conj().T

        return rho
#--------------------------------------------------------------------------------------------------------------------------------------------

#Define a CNOT function with depolarization error:
def DE_fiveQCNOT(control, target, rho, depol_prob=0.01):
    cnot_gates = {
        (0, 1): CNOT_0_1,
        (1, 2): CNOT_1_2,
        (0, 3): CNOT_0_3,
        (1, 3): CNOT_1_3,
        (1, 4): CNOT_1_4,
        (2, 4): CNOT_2_4
    }
    key = (control, target)
    if key in cnot_gates:
        rho = cnot_gates[key] @ rho @ cnot_gates[key].T
    else:
        return 'This CNOT is not available yet'
    rho = depolarizing_error(rho, depol_prob, target)
    return rho


#--------------------------------------------------------------------------------------------------------------------------------------------

def state_to_density_matrix(state):
    """
    Convert a state vector to its corresponding density matrix.

    Parameters:
    state (np.ndarray): The state vector to be converted.

    Returns:
    np.ndarray: The density matrix representation of the state.
    """
    # The outer product of the state with its conjugate to form the density matrix
    return np.outer(state, np.conj(state))

#--------------------------------------------------------------------------------------------------------------------------------------------
def encode_logical_density_matrix(input_density_matrix, depol_prob=0.01):
    initial_density = tensor_product(input_density_matrix, np.outer(ket_0, ket_0), np.outer(ket_0, ket_0))
    encoded = DE_fiveQCNOT(0, 1, initial_density, depol_prob)
    encoded = DE_fiveQCNOT(1, 2, encoded, depol_prob)
    return encoded

#--------------------------------------------------------------------------------------------------------------------------------------------
def rho_add_ancilla_qubits(error_density_matrix, depol_prob=0.01):
    ancilla_density_matrix = tensor_product(error_density_matrix, ket_0 @ ket_0.T, ket_0 @ ket_0.T)
    ancilla_density_matrix = DE_fiveQCNOT(0, 3, ancilla_density_matrix, depol_prob)
    ancilla_density_matrix = DE_fiveQCNOT(1, 3, ancilla_density_matrix, depol_prob)
    ancilla_density_matrix = DE_fiveQCNOT(1, 4, ancilla_density_matrix, depol_prob)
    ancilla_density_matrix = DE_fiveQCNOT(2, 4, ancilla_density_matrix, depol_prob)
    return ancilla_density_matrix

#--------------------------------------------------------------------------------------------------------------------------------------------

def rho_measure_rightmost_2_qubits(rho):
    """
    Measure the rightmost two qubits of a 5-qubit density matrix and return the result as a string.

    Parameters:
    rho (np.ndarray): The density matrix of the 5-qubit system to be measured.

    Returns:
    str: The measurement outcome of the rightmost two qubits as a binary string.
    """
    # Define projection operators for two qubits 
    P_00 = np.kron(ket_0, ket_0).dot(np.kron(ket_0, ket_0).T)
    P_01 = np.kron(ket_0, ket_1).dot(np.kron(ket_0, ket_1).T)
    P_10 = np.kron(ket_1, ket_0).dot(np.kron(ket_1, ket_0).T)
    P_11 = np.kron(ket_1, ket_1).dot(np.kron(ket_1, ket_1).T)

    # Tensor projection operators with identity for the 3 leftmost qubits
    I_3 = np.eye(8)  # Identity for 3 qubits
    P_00_full = np.kron(I_3, P_00)
    P_01_full = np.kron(I_3, P_01)
    P_10_full = np.kron(I_3, P_10)
    P_11_full = np.kron(I_3, P_11)

    # Compute probabilities for each outcome
    p_00 = np.trace(P_00_full @ rho).real
    p_01 = np.trace(P_01_full @ rho).real
    p_10 = np.trace(P_10_full @ rho).real
    p_11 = np.trace(P_11_full @ rho).real

    # print(p_00, p_01, p_10, p_11)
    # Randomly choose a measurement outcome based on the probabilities
    outcomes = ['00', '01', '10', '11']
    probabilities = [p_00, p_01, p_10, p_11]
    chosen_outcome = np.random.choice(outcomes, p=probabilities)


    # Collapse the density matrix onto the chosen outcome
    if chosen_outcome == '00':
        rho_post = (P_00_full @ rho @ P_00_full.T) / p_00
    elif chosen_outcome == '01':
        rho_post = (P_01_full @ rho @ P_01_full.T) / p_01
    elif chosen_outcome == '10':
        rho_post = (P_10_full @ rho @ P_10_full.T) / p_10
    else:  # '11'
        rho_post = (P_11_full @ rho @ P_11_full.T) / p_11


    return chosen_outcome, rho_post    

#--------------------------------------------------------------------------------------------------------------------------------------------

def rho_correct_density_matrix(rho, result, depol_prob=0.01):
    XII = tensor_product(X_gate, Identity, Identity, Identity, Identity)
    IXI = tensor_product(Identity, X_gate, Identity, Identity, Identity)
    IIX = tensor_product(Identity, Identity, X_gate, Identity, Identity)
    if result == '00':
        return rho
    elif result == '01':
        rho = IIX @ rho @ IIX.T
        return depolarizing_error(rho, depol_prob, 2)
    elif result == '10':
        rho = XII @ rho @ XII.T
        return depolarizing_error(rho, depol_prob, 0)
    elif result == '11':
        rho = IXI @ rho @ IXI.T
        return depolarizing_error(rho, depol_prob, 1)
    else:
        raise ValueError("Invalid measurement result. Expected one of ['00', '01', '10', '11'].")


#--------------------------------------------------------------------------------------------------------------------------------------------

def fidelity(rho, sigma):
    """
    Compute the fidelity between two density matrices, which quantifies the closeness of two quantum states.

    Parameters:
    rho (np.ndarray): The first density matrix, representing a quantum state.
    sigma (np.ndarray): The second density matrix, representing a quantum state.

    Returns:
    float: The fidelity, ranging from 0 to 1, where 1 indicates identical states.
    """
    # Compute the square root of the first density matrix
    sqrt_rho = sqrtm(rho)
    # Form the product of sqrt(rho), sigma, and sqrt(rho)
    product = np.dot(sqrt_rho, np.dot(sigma, sqrt_rho))
    # Compute the square root of the product matrix
    sqrt_product = sqrtm(product)
    # The fidelity is the square of the trace of sqrt_product
    return np.real((np.trace(sqrt_product))**2) # Taking the real part to avoid any small imaginary parts due to numerical errors

#--------------------------------------------------------------------------------------------------------------------------------------------
