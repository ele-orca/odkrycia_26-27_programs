#code by chatgpt so I dont know what is the function dispresion(k, f), I just know it works for an infinite ferromagnet magntized prepedicularlly to its surface
#if I learn anything I will write about it here.
from scipy.optimize import brentq
import numpy as np

#constants
Aex = 18.5e-12
Msat = 1.1e6
mu0 = 4 * np.pi * 1e-7
gamma = 28e9 * 2 * np.pi
Ku = 8.3e5

f = 10e9
thickness = 1e-9

D = 2 * Aex / Msat
Bi = 2 * Ku / Msat - mu0 * Msat


def F(a):
    x = a * thickness
    return 1 - (1 - np.exp(-x)) / x


def dispersion(k, f):
    B = Bi + D*k**2

    return (
        gamma*np.sqrt(
            B*(B + mu0*Msat*F(k))
        ) - 2*np.pi*f
    )


k = brentq(lambda k: dispersion(k,f), 1e5, 1e9)
wavelength = 2 * np.pi / k

print(wavelength)
