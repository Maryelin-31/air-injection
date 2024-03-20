# -*- coding: utf-8 -*-

from zml import *
from zmlx.fluid import *
from zml import Interp2, TherFlowConfig, data_version

'gas mixture'
"gas mixture density, viscosity the same C2h6(ethane)"
from Utility.gas_density.gas_mixture_density import GAS_den
from zmlx.fluid.conf.gas_viscosity.c2h6 import gas_vis_c2h6


import os
import math



"Gas Mixture"

def Gas_mixture(tmin=280, tmax=2000, pmin=1.0e6, pmax=40.0e6):
    def gas_den(P, T):
        density = GAS_den(P, T)
        return density

    def get_density(P, T):
        return gas_den(P, T)

    def create_density():
        den = Interp2()
        den.create(pmin, 1e6, pmax, tmin, 10, tmax, get_density)
        return den

    def gas_vis(P, T):
        viscosity = gas_vis_c2h6(P, T)
        return viscosity

    def get_viscosity(P, T):
        return gas_vis(P, T)

    def create_viscosity():
        vis = Interp2()
        vis.create(pmin, 1e6, pmax, tmin, 10, tmax, get_viscosity)
        return vis
    specific_heat = 1800

    return TherFlowConfig.FluProperty(den=create_density(), vis=create_viscosity(), specific_heat=specific_heat)


