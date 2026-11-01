import numpy as np

"""this function can be used to choose the potential(pot)
pot: potential of free energy
para: parameter alpha, chi or H depending in the potential"""
def df_dc(c, potential, parameter=1.0):

    if potential == "asymmetric":
        return c**3-c+parameter

    elif potential == "flory_huggins":
        c_safe = np.clip(c, 1e-8, 1 - 1e-8)
        return np.log(c_safe / (1 - c_safe)) + parameter * (1 - 2 * c_safe)
#np.clip() prevents taking the logarithm of zero or a negative value
    elif potential == "triple_well":
        return 2 * c * (c ** 2 - 1) * (3 * c ** 2 - 1)
    elif potential == "external_field":
        return c ** 3 - c - parameter

    elif potential == "periodic":
        return np.sin(c)

    else:
        raise ValueError("Unknown potential: " + potential)

#test
"""if __name__ == "__main__":
    c = np.array([[-1.0, 0.0, 1.0]])

    print("Asymmetric:", df_dc(c, "asymmetric", 1.0))
    print("Triple well:", df_dc(c, "triple_well"))
    print("External field:", df_dc(c, "external_field", 1.0))
    print("Periodic:", df_dc(c, "periodic"))
"""