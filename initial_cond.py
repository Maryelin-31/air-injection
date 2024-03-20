# -*- coding: utf-8 -*-

# from Utility.Chang7.get_path import *
# from zmlx.alg.clamp import clamp
# from saturation_data import *
from zml import *


def create_initial():
    """
    create initial field
    """

    def get_initial_t(x, y, z):
        """
        the initial temperature
        """
        return 338.0 + 22.15 - 0.0443 * y

    def get_initial_p(x, y, z):
        """
        the initial pressure
        """
        return 10e6 - 1e4 * y
    
    def get_perm(x, y, z):
        """
        the initial permeability
        """
        # t = Tensor3()
        # t.xx = 1.0e-14
        # t.yy = 1.0e-15
        # t.zz = 1.0e-15
        return 1.0e-15

    def get_initial_s(x, y, z):
        """
        the initial saturation ()
        """

        return (0, 0, 0, 0, 0, 0), 0, 0, 0, (1.0, 0)
    
    
    def get_fai(x, y, z):
        """
        porosity
        """
        return 0.43

    def get_denc(x, y, z):
        """
        density * heat capacity
        """
        return 2600 * 1000

    def get_heat_cond(x, y, z):
        return 1.0

    return {'porosity': get_fai, 'pore_modulus': 50e6, 'p': get_initial_p,
            'temperature': get_initial_t,
            'denc': get_denc, 's': get_initial_s,
            'perm': get_perm, 'heat_cond': get_heat_cond, 'dist': 0.01}


if __name__ == '__main__':
    print(create_initial())
