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
from zmlx.utility.PressureController import PressureController
from zmlx.utility.SeepageCellMonitor import SeepageCellMonitor

# from zml_hyd import *
from matplotlib import pyplot as pyplot
from matplotlib import ticker, cm
from matplotlib.ticker import LinearLocator
from matplotlib.cm import ScalarMappable
import os
import numpy as np




# Create seepage mesh.

def create_mesh():
    mesh = SeepageMesh.create_cube(np.linspace(0, 50, 101), np.linspace(0, 50, 26), (-1.5, 1.5))
    return mesh

def create_reaction():
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
    
    config.reactions.append(combustion.create(left=[(config.iwat, 0.4), ((config.iGAS, 5), 0.60)], 
                                              right=[((config.isol, 0), 0.65), ((config.iGAS, 3), 0.35)],
                            temp = 623.15, 
                            heat = 783e6, 
                            rate = 1.0e-6),)
    
    config.reactions.append(combustion.create(left=[((config.iGAS, 4), 1.0)], 
                                              right=[((config.iGAS, 5), 0.8), ((config.iGAS, 3), 0.2)],
                            temp = 623.15, 
                            heat = 0.305e6, 
                            rate = 1.0e-6),)
    
    config.reactions.append(combustion.create(left=[(config.iwat, 0.25), ((config.iGAS, 4), 0.75)], 
                                              right=[((config.isol, 1), 0.6), ((config.iGAS, 3), 0.4)],
                            temp = 623.15, 
                            heat = 1.0e6, 
                            rate = 1.0e-6),)
 
    # h2o and steam
    config.reactions.append(
        vapor_react.create(
            vap=(config.iGAS, 2),
            wat=config.iwat,
            fa_t=config.flu_keys['temperature'],
            fa_c=config.flu_keys['specific_heat']))
 
    # The decomposition of Kerogen.
    config.reactions.append(
        decomposition.create(left=(config.isol, 0), right=[(config.iHO, 0.663), 
                                                            (config.iLO, 0.076),
                                                            ((config.iGAS, 0), 0.046),
                                                            ((config.isol, 1), 0.215),
                                                            ],
                             temp=563.15, heat=161100.0, 
                             rate=4.81e-6,
                             fa_t=config.flu_keys['temperature'],
                             fa_c=config.flu_keys['specific_heat']))
    # The decomposition of Heavy oil
    config.reactions.append(
        decomposition.create(left=config.iHO, right=[(config.iLO, 0.438),   
                                                         ((config.iGAS, 1), 0.217),
                                                         ((config.isol, 1), 0.345),
                                                         ],
                             temp=623.15, heat=219328.0, 
                             rate=2.71e-7,
                             fa_t=config.flu_keys['temperature'],
                             fa_c=config.flu_keys['specific_heat']))
    return config


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
        return 15.0e6 + 5e6 - 1e4 * y

    def get_perm(x, y, z):
        """
        the initial permeability
        """
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
    
    return {'porosity': get_fai, 'pore_modulus': 100e6, 'p': get_initial_p,
            'temperature': get_initial_t,
            'denc': get_denc, 's': get_initial_s,
            'perm': get_perm, 'heat_cond': get_heat_cond, 'dist': 0.01}


"Create Mesh"
mesh = create_mesh()
# print(mesh)
    
config = create_reaction()

"Initial Conditions Reservoir"
ini = create_initial()

"create model"
model = config.create(mesh, **ini)  # create
model.set(gravity=(0, -10, 0))


"relative permeability"


def stone_model_I(swir, sorg, sorw, sgc, krwro, kroiw, krgro, nw, nsorw, ng, nog):
    assert swir < 1
    #oil-water system and gas-oil system Corey two phases model
    #variables
    sw = np.linspace(swir, 1 - sorw, 20, endpoint=True)
    sg = np.linspace(sgc, 1 - sorg, 20, endpoint=True)
    so = 1 - sg
    #Models Corey, 1954
    krw = krwro * ((sw - swir) / (1 - sorw - swir))**nw

    krow = kroiw * ((1 - sw - sorw) / (1 - sorw - swir))**nsorw

    krg = krgro * ((sg - sgc) / (1 - sgc - sorg - swir))**ng
    krg[krg >=1] = 1

    krog = kroiw * ((1 - sg - sorg - swir) / (1 - sgc - sorg - swir))**nog

    #Stone Model I normalized by Aziz and Settari, 1979
    #swc = swir
    #Fayers and Mattews 1984
    a = 1 - (sg / (1 - swir - sorg))
    som= (a * sorw) + ((1 - a) * sorg)
    s_o = np.abs(so - som) / (1 - swir - som)  # so>= som
    s_w = np.abs(sw - swir) / (1 - swir - som)  # sw >= swir
    s_g = (sg) / (1 - swir - som)
    s_o[s_o >= 1.0] = 1 - swir
    s_w[s_w >= 1.0] = 1 - sorw
    s_g[s_g >= 1.0] = 1 - sorg
    kro0 = kroiw
    kro = (s_o / kro0) * (krow / (1 - s_w)) * (krog / (1 - s_g))
    kro[kro >= 1] = 1
    return sw, krw, sg, krg, so, kro
sw, krw, sg, krg, so, kro = stone_model_I(swir=0.1, sorg=0.1, sorw=0.1, sgc=0.1, 
                                          krwro=0.9, kroiw=1, krgro=0.9, 
                                          nw=2, nsorw=2, ng=2, nog=2)

# Set relative permeability.
model.set_kr(config.iGAS, sg, krg)
model.set_kr(config.iwat, sw, krw)
model.set_kr(config.iHO, so, kro)


"Production Well"

pos = (25, 25, 1)
cell = model.get_nearest_cell(pos=pos)
x, y, z = cell.pos
virtual_cell = config.add_cell(model, vol=1.0e6 , pos=pos, porosity=0.1, pore_modulus=20e6,
                               temperature=ini['temperature'](*pos), p=3.5e6,
                               s=ini['s'](*pos))

face = config.add_face(model, virtual_cell, model.get_cell(cell.index),
                       heat_cond=0, perm=ini['perm'](*pos), area=3.0, length=1.0, )

lenght = face.get_attr(config.face_keys['length'])
cond = (2 * np.pi * z *  ini['perm'](*pos)) / (np.log((lenght)**2) - np.log(3.0)) #radial flow
face.set_attr(config.face_keys['g0'], cond)

prectrl = PressureController(virtual_cell, t=[0, 1e10], p=[3.5e6, 3.5e6])
monitor = SeepageCellMonitor(get_t=lambda: config.get_time(model),
                              cell=(virtual_cell, prectrl))


"Injection Configuration"
rate_inj = 0.002
cell_ids = set()

#Injection Temperature
temp_inje = 900
z = 0
y = 25.0
x = 12.25
# for x in (12.5, 30):
cell = model.get_nearest_cell(pos=(x, y, z))
cell_ids.add(cell.index)

# Add_injector

for cell_id in cell_ids:
    fa_t = 0
    fa_c = 1
    cell1 = model.get_cell(cell_id)
    flu = cell1.get_fluid(0).get_component(3)  # Air
    flu.set_attr(fa_t, 700)
    flu.set_attr(fa_c, 1000)
    model.add_injector(cell=cell1, fluid_id=[0, 3], flu=flu,
                        pos=cell1.pos, radi=1.0, opers=[(0, rate_inj)])

"Time step Strategy"
config.set_dv_relative(model, 0.5)  # The ratio of the distance traveled by each time step to the grid size
config.set_dt(model, 0.01)  # initial value for time step
config.set_dt_max(model, 24 * 3600 * 7)  # Maximum value of time step <one week> # Maximum value of time step <one week>
config.set_dt_min(model, 3600)  # Minimum step size is 1 hour


"Save Results"


    
def cell_mass(cell):
    fluid = []
    total = []
    for i in range(cell.fluid_number):
        total.append(cell.get_fluid(i).mass)
        if cell.get_fluid(i).component_number == 0:
            fluid.append(cell.get_fluid(i).mass)
        else:
            for j in range(cell.get_fluid(i).component_number):
                fluid.append(cell.get_fluid(i).get_component(j).mass)
    saturation = [i / sum(total) for i in fluid]
    return saturation

def save_wt(path):
    name = os.path.basename(__file__)
    result_folder = os.path.join(os.getcwd(), f'data_{name}', f'Results_rate', 'wt_cells')
    
    if os.path.exists(f'result_folder'):
        import shutil
        shutil.rmtree(f'result_folder')            
    os.makedirs(result_folder, exist_ok=True)
    
    SavePath = os.path.join(result_folder, path)
    with open(SavePath, 'w') as file:
        for cell in model.cells:
            x, y, z = cell.pos
            satu = cell_mass(cell)
            satu_str = ' '.join(str(i) for i in satu)
            file.write(f'{x} {y} {z} {satu_str}\n')


def save_mass(path):
    name = os.path.basename(__file__)
    result_folder = os.path.join(os.getcwd(), f'data_{name}', f'Results_rate', 'Results_cells')
    
    if os.path.exists(f'result_folder'):
        import shutil
        shutil.rmtree(f'result_folder')            
    os.makedirs(result_folder, exist_ok=True)
    
    SavePath = os.path.join(result_folder, path)
    with open(SavePath, 'w') as file:
        for cell in model.cells:
            x, y, z = cell.pos
            temp = cell.get_attr(config.cell_keys['temperature'])
            pres = cell.get_attr(config.cell_keys['pre'])
            file.write(f'{x} {y} {z} '
                       f'{cell.get_fluid(0).get_component(0).mass} {cell.get_fluid(0).get_component(1).mass} {cell.get_fluid(0).get_component(2).mass} '
                       f'{cell.get_fluid(0).get_component(3).mass} {cell.get_fluid(0).get_component(4).mass} {cell.get_fluid(0).get_component(5).mass} '
                       f'{cell.get_fluid(1).mass} '
                       f'{cell.get_fluid(2).mass} '
                       f'{cell.get_fluid(3).mass} '
                       f'{cell.get_fluid(4).get_component(0).mass} {cell.get_fluid(4).get_component(1).mass} '
                       f'{temp} {pres}\n')

def fluid_sat(cell):
    fluid = []
    for i in range(cell.fluid_number):
        if cell.get_fluid(i).component_number == 0:
            fluid.append(cell.get_fluid(i).vol_fraction)
        else:
            for j in range(cell.get_fluid(i).component_number):
                fluid.append(cell.get_fluid(i).vol_fraction)
    return fluid

def save_sat(path):
    name = os.path.basename(__file__)
    result_folder = os.path.join(os.getcwd(), f'data_{name}', f'Results_rate', 'Saturation_cells')
    
    if os.path.exists(f'result_folder'):
        import shutil
        shutil.rmtree(f'result_folder')            
    os.makedirs(result_folder, exist_ok=True)
    
    SavePath = os.path.join(result_folder, path)
    with open(SavePath, 'w') as file:
        for cell in model.cells:
            x, y, z = cell.pos
            temp = cell.get_attr(config.cell_keys['temperature'])
            pres = cell.get_attr(config.cell_keys['pre'])
            satu = fluid_sat(cell)
            satu_str = ' '.join(str(i) for i in satu)
            file.write(f'{x} {y} {z} {satu_str} {temp} {pres}\n')

def run_2():
    name = os.path.basename(__file__)
    
    data_path = os.path.join(os.getcwd(), f'data_{name}')
    resu_path = os.path.join(data_path, f'Results_rate')
    os.makedirs(os.path.join(os.getcwd(), data_path, f'Results_rate'), exist_ok=True)

    # folder = os.getcwd()
    folder = os.path.join(os.getcwd(), data_path, f'Results_rate')

    
    solver = ConjugateGradientSolver()
    solver.set_tolerance(1.0e-8)
    for step in range(1000000):
        config.iterate(model, solver=solver)
        time = config.get_time(model)
        
        if config.get_time(model) > 3600 * 24 * 365 * 20:
            print(f'{step}, Finish')
            break
        
        path = f'time_{time}.txt'
        if step % 100 == 0:
            save_wt(path)
            save_sat(path)
            save_mass(path)
            monitor.update(dt=config.get_dt(model))
            monitor.save(os.path.join(folder, f'prod_{name}.txt'))
            print(f'{step} {time / (3600 * 24 * 365)}')
    
# run_2()    
