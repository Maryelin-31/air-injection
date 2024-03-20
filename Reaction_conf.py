# -*- coding: utf-8 -*-
from zml import *
from zmlx.fluid import *
from zmlx.fluid.conf import *
from zmlx.config.TherFlowConfig import TherFlowConfig
from zmlx.fluid.kerogen import create as create_kerogen
from zmlx.fluid.h2o_gas import create as create_h2o_gas
from zmlx.fluid.h2o import create as create_h2o
from zmlx.fluid.char import create as create_char
from zmlx.fluid.ch4_lyx import create as create_methane_gas
from zmlx.fluid.conf.CO2_gas import create_flu as create_co2
from zmlx.fluid.conf.gas_mixture import Gas_mixture
from zmlx.fluid.conf.O2_gas import create_flu as create_o2
from zmlx.fluid.conf.CO_gas import create_flu as create_co


from zmlx.kr.create_krf import create_krf
from zmlx.react import decomposition
from zmlx.react import combustion
from zmlx.react import vapor as vapor_react


def create():
    """

    """
    config = TherFlowConfig()
    
    # Gas Phase (gas mixture, methane, steam_water())
    config.iGAS = config.add_fluid([Gas_mixture(), create_methane_gas(), create_h2o_gas(), create_o2(), create_co2(), create_co()])
    
    # Water
    config.iwat = config.add_fluid(create_h2o())
    
    #Light Oil
    config.iLO = config.add_fluid(create_light_oil_liq())
    
    # Heavy Oil
    config.iHO = config.add_fluid(create_heavy_oil_liq())

    # Solid phase Kerogen
    config.isol = config.add_fluid([create_kerogen(), create_char()])
    
    # Kerogen Combustion
    
    config.reactions.append(combustion.create(left=[(config.iwat, 0.2), ((config.iGAS, 4), 0.5), ((config.iGAS, 5), 0.3)], 
                                              right=[((config.isol, 0), 0.6), ((config.iGAS, 3), 0.4)],
                            temp = 623.15, 
                            heat = 1.640e6, 
                            rate = 1.0e-2),)

    # h2o and steam
    config.reactions.append(
        vapor_react.create(
            vap=(config.iGAS, 2),
            wat=config.iwat,
            fa_t=config.flu_keys['temperature'],
            fa_c=config.flu_keys['specific_heat']))

    # The decomposition of Kerogen.
    config.reactions.append(
        decomposition.create(left=(config.isol, 0), right=[(config.iHO, 0.6), 
                                                            (config.iLO, 0.1),
                                                            (config.iwat, 0.05),
                                                            ((config.iGAS, 0), 0.05),
                                                            ((config.isol, 1), 0.2),
                                                            ],
                             temp=563.15, heat=161100.0, 
                             rate=4.81e-6,
                             fa_t=config.flu_keys['temperature'],
                             fa_c=config.flu_keys['specific_heat']))

    # The decomposition of Heavy oil
    config.reactions.append(
        decomposition.create(left=config.iHO, right=[(config.iLO, 0.5),   
                                                         ((config.iGAS, 0), 0.2),
                                                         ((config.isol, 1), 0.3),
                                                         ],
                             temp=623.15, heat=219328.0, 
                             rate=2.71e-7,
                             fa_t=config.flu_keys['temperature'],
                             fa_c=config.flu_keys['specific_heat']))
    

    return config


if __name__ == '__main__':
    c = create()
