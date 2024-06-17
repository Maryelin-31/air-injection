from zml import *
from zmlx.config import seepage
from zmlx.seepage_mesh.cube import create_xz

from zmlx.fluid.ch4 import create as create_ch4
from zmlx.fluid.h2o_gas import create as create_steam
from zmlx.fluid.conf.O2_gas import create as create_o2
from zmlx.fluid.conf.CO2_gas import create as create_co2
from zmlx.fluid.conf.CO_gas import create as create_co
from zmlx.fluid.h2o import create as create_h2o
from zmlx.fluid.conf.c11h24_liq import create as create_lo
from zmlx.fluid.conf.c22h46_liq import create as create_ho
from zmlx.fluid.kerogen import create as create_kerogen
from zmlx.fluid.char import create as create_char

from zmlx.react import decomposition
from zmlx.react import combustion
from zmlx.react import vapor
from zmlx.react.inh import add_inh

from zmlx.kr.create_krf import create_krf
from zml import Interp1, create_dict
from zmlx.utility.PressureController import PressureController
from zmlx.utility.SeepageCellMonitor import SeepageCellMonitor


def create_mesh():
    mesh = create_xz(x_min=0, dx=0.5, x_max=100,
                     y_min=-1, y_max=0,
                     z_min=-100, dz=0.5, z_max=0,)
    return mesh

def create_fludefs():
    """
    0. gas = [methane, steam, oxigen, co, co2]
    1. h20
    2. lo = light oil
    3. ho = heavy oil
    4. sol= kerogen, char
    
    """
    
    gas = Seepage.FluDef(name='gas')
    gas.add_component(create_ch4(name='ch4'))
    gas.add_component(create_steam(name='steam'))
    gas.add_component(create_o2(name='o2'))
    gas.add_component(create_co(name='co'))
    gas.add_component(create_co2(name='co2'))
    
    h2o = create_h2o(name='h2o')
    lo = create_lo(name='lo')
    ho = create_ho(name='ho')
    
    sol = Seepage.FluDef(name='sol')
    sol.add_component(create_kerogen(name='kero'))
    sol.add_component(create_char(name='coke'))
    
    return [gas, h2o, lo, ho, sol]

def create_reactions(temp_max=None):
    
    """
    1. Steam phase transitions
    2. Combustion = 10.1016/j.egyr.2022.08.174: 
        Table 1, number 1, 9, 10
        Expresses in mass fractions 
        0.65 kero + 0.35 O2 --> 0.4H2O + 0.6CO + ENERGY
        0.8CO + 0.2 O2 --> CO2 + Energy
        0.6coke + 0.4 O2 --> 0.25H2O + 0.4CO2 + Energy
        
    3. Kerogen Decomposition:
        a simplification of 10.2172/10169154
        Table 1
    """
    
    result = []
    
    r = vapor.create(vap='steam', wat='h2o', temp_max=temp_max)
    result.append(r)
    
    r = combustion.create(left =[('h2o', 0.4), ('co', 0.6)],
                          right=[('kero', 0.65), ('o2', 0.35)],
                          temp = 623.15, heat=783e6, rate=1)
    result.append(r)
    r = combustion.create(left =[('co2', 1)],
                          right=[('co', 0.8), ('o2', 0.2)],
                          temp = 623.15, heat=0.305e6, rate=1)
    result.append(r)
    r = combustion.create(left =[('h2o', 0.25), ('co', 0.75)],
                          right=[('coke', 0.6), ('o2', 0.4)],
                          temp = 623.15, heat=1.0e6, rate=1)
    result.append(r)
    
    r = decomposition.create(left ='kero',
                             right=[('ho', 0.663), ('lo', 0.076), ('ch4', 0.046), ('coke', 0.215)], 
                             temp=563.15, heat=161100.0, rate=4.81e-6)
    result.append(r)
    r = decomposition.create(left ='ho',
                             right=[('lo', 0.438), ('ch4', 0.217), ('coke', 0.345)], 
                             temp=623.15, heat=219328.0, rate=2.71e-7)
    # When the solid occupies 80% of the total pores, increase the cracking temperature to limit further decomposition (avoid all pores being occupied by solids)
    add_inh(r, sol='sol', liq=None,
            c=[0, 0.8, 1.0],
            t=[0, 0, 1.0e4])
    result.append(r)
    return result

def create_initial():
    """
    create initial field
    """

    def get_initial_t(x, y, z):
        """
        the initial temperature
        """
        return 338.0 + 22.15 - 0.0443 * z

    def get_initial_p(x, y, z):
        """
        the initial pressure
        """
        return 15.0e6 + 5e6 - 1e4 * z

    def get_perm(x, y, z):
        """
        the initial permeability
        """
        return 1.0e-15

    def get_initial_s(x, y, z):
        """
        the initial saturation ()
        """

        return (0.08, 0, 0, 0, 0), 0.04, 0.08, 0.2, (0.6, 0)
    
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
    
    return {'porosity': get_fai, 'pore_modulus': 100e6, 'p': get_initial_p,
            'temperature': get_initial_t,
            'denc': get_denc, 's': get_initial_s,
            'perm': get_perm, 'heat_cond': get_heat_cond, 'dist': 0.01}


x, y = create_krf(faic=0.02, n=3.0, k_max=100,
                  s_max=2.0, count=500)

gr = Interp1(x=x, y=y)

kw = create_dict(fludefs=create_fludefs(),
                 reactions=create_reactions(),
                 gr=gr,
                 has_solid=False,)


kw.update(**create_initial())

gravity = [0, 0, -10]


kw.update(create_dict(dt_max=3600.0 * 24.0 * 10.0,
                      gravity=gravity, ))

model = seepage.create(mesh=create_mesh(), **kw)
