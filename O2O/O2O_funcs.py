# Math, numerics, and graphing
import numpy as np
import scipy as sp
import scipy.integrate as integrate
from scipy.optimize import minimize
from scipy.special import factorial
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from scipy.integrate import trapz
from scipy.special import eval_hermite
from scipy.optimize import minimize_scalar
import strawberryfields as sf

def fock_to_position(n, q):
    '''Returns the Fock position wavefunctions defined over q space.

    Args:
        n (int): Fock index
        q (array): position values

    Returns:
        position_wf_fock (array): nth Fock state wavefunction
    '''

    position_wf_fock = ((1/np.sqrt((2**n)*factorial(n)))*(np.pi**-0.25)
            *np.exp(-(q**2)/2)*eval_hermite(n,q) ) #Employs the scipy.special eval_hermite function

    return position_wf_fock

def fock_coefficients_to_position(coeffs, q):
    '''Returns the position wavefunction for a superposition of Fock states
    defined by the given Fock coefficients.

    Args:
        coeffs (array): array of Fock coefficients
        q (array): position values over which to compute the wavefunction

    Returns:
        wavefunction (array): resulting wavefunction in the position basis
    '''
    wavefunction = np.zeros_like(q, dtype=complex)  # Initialize the wavefunction

    # Loop over all Fock coefficients and sum the weighted Fock states
    for n, coeff in enumerate(coeffs):
        fock_wf = fock_to_position(n, q)  # Get the nth Fock state in position basis
        wavefunction += coeff * fock_wf   # Add the weighted Fock state to the total wavefunction

    return wavefunction

def qunaught_fock_coeff(n, combs=100, alpha=np.sqrt(2*np.pi), aspect_ratio=1):
    '''Returns the nth Fock coefficient for the ideal qunaught GKP state.

    Args:
        n (int): index of Fock state
        combs (int): number of delta spikes in the ideal state used to 
                        calculate inner product with Fock state
        alpha (float): Qunaught lattice constant (default is sqrt(2*pi))
        aspect_ratio (float): the aspect ratio lambda that scales position and momentum lattice spacing

    Returns:
        coeff (complex): inner product between ideal qunaught GKP and nth 
            Fock states.
    '''
    # adjust alpha by aspect_ratio in the q direction
    alpha_q = alpha * aspect_ratio
    # q values from which to sample the nth Fock state
    samples = (np.arange(-combs/2, 1+combs/2))*alpha_q
    # sum of sampled values to yield inner product
    coeff = np.sum(fock_to_position(n, samples))
    return coeff


def qunaught_fock(eps, q, n_max=100, norm=True, aspect_ratio=1): 
    '''Returns the epsilon-qunaught q wavefunction and coefficients in the 
        Fock basis.

    Args:
        eps (float): squeezing parameter (epsilon) for the qunaught state
        q (array): array of q values for q wavefunction
        n_max (int): Fock cutoff (maximum number of Fock states considered)
        norm (Boolean): whether to return a normalized state
        aspect_ratio (float): the aspect ratio lambda that scales position and momentum lattice spacing

    Returns:
        qunaught (array): q wavefunction for qunaught state
        coeffs (array): Fock coefficients up to cutoff.
    '''
    alpha = np.sqrt(2*np.pi)  # Lattice constant for qunaught is sqrt(2*pi)

    qunaught = 0
    coeffs = np.zeros(n_max+1, dtype=complex)
    # initialize normalization constant
    N = 0
    for i in range(n_max+1):
        # calculate nth coefficient for the qunaught state, with aspect_ratio
        coeff = qunaught_fock_coeff(i, alpha=alpha, aspect_ratio=aspect_ratio) * np.exp(-i * eps)
        coeffs[i] = coeff
        qunaught += coeff * fock_to_position(i, q)
        if norm:
            N += np.absolute(coeff)**2
    
    # Check if normalization constant is non-zero before dividing
    if norm and N > 0:
        qunaught = qunaught / np.sqrt(N)
        coeffs = coeffs / np.sqrt(N)
    elif N == 0:
        print("Warning: Normalization constant is zero, skipping normalization.")
        
    return qunaught, coeffs

def agn_sample(variance):
    """
    Samples a single displacement from an Additive Gaussian Noise (AGN) channel.
    """
    if variance < 0:
        raise ValueError("Variance must be non-negative.")
    
    # Generate a single Gaussian sample with mean = 0 and variance = sigma^2
    sample = np.random.normal(loc=0.0, scale=np.sqrt(variance))
    return sample

# Define the modular reduction function using sf.math
def reduce_modulo(value, modulus=np.sqrt(2 * np.pi)):
    """
    Reduce the input value modulo the specified modulus.
    """
    nearest_multiple = sf.math.floor(value / modulus + 0.5)  # Closest integer to value/modulus
    reduced_value = value - nearest_multiple * modulus
    return reduced_value

def decode_q(measurement, modulus=np.sqrt(2 * np.pi)):
    """
    Compute the redisplacement value for the position quadrature.
    """
    reduced_value = reduce_modulo(measurement, modulus) 
    z1 = 0 ##How to obtain Z1?
    redisplacement = -0.5 * (reduced_value + z1)
    return redisplacement

def decode_p(measurement, modulus=np.sqrt(2 * np.pi)):
    """
    Compute the redisplacement value for the momentum quadrature.
    """
    reduced_value = reduce_modulo(measurement, modulus)
    redisplacement = reduced_value
    return redisplacement