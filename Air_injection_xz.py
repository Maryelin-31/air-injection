import numpy as np
from matplotlib import pyplot as pyplot
from matplotlib.cm import ScalarMappable
import os

from zml import *
from zmlx.config import seepage
from zmlx.seepage_mesh.cube import create_xz
from zml import Dfn2
from zmlx.plt.show_dfn2 import show_dfn2
from zmlx.geometry.seg_point_distance import seg_point_distance
from zmlx.geometry.point_distance import point_distance

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
from zmlx.utility.SaveManager import SaveManager

def create_mesh():
    mesh = create_xz(x_min=0, dx=0.5, x_max=100,
                     y_min=-1, y_max=0,
                     z_min=0, dz=0.5, z_max=100,)
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
    2. Combustion = doi:10.1016/j.egyr.2022.08.174: 
        Table 1, number 1, 9, 10
        Expresses in mass fractions 
        0.65 kero + 0.35 O2 --> 0.4H2O + 0.6CO + ENERGY
        0.8CO + 0.2 O2 --> CO2 + Energy
        0.6coke + 0.4 O2 --> 0.25H2O + 0.4CO2 + Energy
        
    3. Kerogen Decomposition:
        a simplification of doi:10.2172/10169154
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
        if 25 <= z <= 75:
            return 1.0e-15  
        else:
            return 0
        
    def get_initial_s(x, y, z):
        """
        the initial saturation ()
        """
        if 25 <= z <= 75:
            return (0.08, 0, 0, 0, 0), 0.04, 0.08, 0.2, (0.6, 0)
        else:
            return (1.0, 0, 0, 0, 0), 0, 0, 0, (0, 0)

    def get_fai(x, y, z):
        """
        porosity
        """
        if 25 <= z <= 75:
            return 0.43
        else:
            return 0.01
        

    def get_denc(x, y, z):
        """
        density * heat capacity
        """
        if 25 <= z <= 75:
            return 2600 * 1000
        else:
            return 2600 * 3000

    def get_heat_cond(x, y, z):
        if 25 <= z <= 75:
            return 1.0
        else:
            return 0.01
    
    
    return {'porosity': get_fai, 'pore_modulus': 100e6, 'p': get_initial_p,
            'temperature': get_initial_t,
            'denc': get_denc, 's': get_initial_s,
            'perm': get_perm, 'heat_cond': get_heat_cond, 'dist': 0.01}

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

"Define Model"
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

"Relative Permeability"

model.set_kr(index=0, saturation=sg, kr=krg)
model.set_kr(index=1, saturation=sw, kr=krw)
model.set_kr(index=2, saturation=so, kr=kro)
model.set_kr(index=3, saturation=so, kr=kro)


"Fractures"
"""
fl_min, fl_max: minimum and maximum lengths of natural fractures
p21: Density of natural fractures
"""

dfn = Dfn2()
dfn.range = [0.25, 25.0, 100, 75]

fl_min, fl_ma = 10, 15
p21 = 0.2
angles=None
dfn.add_frac(angles=np.linspace(0.0, 3.1415 * 2, 50) if angles is None else angles,
              lengths=np.linspace(fl_min, fl_ma, 50), p21=p21)
fractures = dfn.get_fractures()
# show_dfn2(fractures, caption='裂缝')
# add cracks to the model
for x0, z0, x1, z1 in dfn.get_fractures():
    cell_beg = model.get_nearest_cell(pos=[x0, 0, z0])
    cell_end = model.get_nearest_cell(pos=[x1, 0, z1])
    
    def get_dist(cell_pos):
        return seg_point_distance([[x0, z0], [x1, z1]], cell_pos[0: 2]) + point_distance(cell_pos, cell_end.pos)
    count = 0
    while cell_beg.index != cell_end.index:
        dist = [get_dist(c.pos) for c in cell_beg.cells]
        idx = 0
        for i in range(1, len(dist)):
            if dist[i] < dist[idx]:
                idx = i
        cell = cell_beg.get_cell(idx)
        face = model.add_face(cell_beg, cell)
        seepage.set_face(face=face, perm=1.0e-14)
        count += 1
        cell_beg = cell


"attribute"
ca = seepage.cell_keys(model)

"Boundaries"

z1, z2 = model.get_pos_range(2)
cells_z1 = model.get_cells_in_range(zr=(z1 - 25, z1 + 25))
cells_z2 = model.get_cells_in_range(zr=(z2 - 25, z2 + 25))
tot_cells= cells_z1 + cells_z2
for cell in tot_cells:
    c = model.get_cell(cell.index)
    mc= c.set_attr(ca.mc, 1.0e10)
    ct= c.set_attr(ca.temperature, 300) 
        
"Injection"
rate_inj = 1.0
pos_inj = (25, 1.0e3, 50)
id_inj  = model.get_nearest_cell(pos=(25, 1.0e3, 50)).index
fa_t = 1
fa_c = 2
cell_inj = model.get_cell(id_inj)
flu = cell_inj.get_fluid(0).get_component(3)
flu.set_attr(fa_t, 700)
flu.set_attr(fa_c, 1000)
model.add_injector(cell=cell_inj, fluid_id=[0, 3], flu=flu,
                    pos=cell_inj.pos, radi=1.0, opers=[(0, rate_inj)])

"Production"

pos_prd = (75, 1e3, 50)
p_prod  = 5e6
id_prod = model.get_nearest_cell(pos_prd).index
virtual_cell = seepage.add_cell(model, pos=pos_prd, porosity=1.0e5, pore_modulus=100e6, vol=1.0,
                                temperature=350,
                                p=p_prod, 
                                s=((1.0, 0, 0, 0, 0), 0, 0, 0, (0, 0)))
seepage.add_face(model, virtual_cell, model.get_cell(id_prod),
                 heat_cond=0, perm=1.0e-14,
                 area=0.0, #the area 0 indicates the well is closed
                 length=1.0)
pre_ctrl = PressureController(virtual_cell, t=[0, 1e10], p=[p_prod, p_prod])
monitor  = SeepageCellMonitor(get_t=lambda: seepage.get_time(model), cell=(virtual_cell, pre_ctrl))

"save"

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
            temp = cell.get_attr(seepage.cell_keys(model).temperature)
            pres = cell.pre
            file.write(f'{x} {y} {z} '
                        f'{cell.get_fluid(0).get_component(0).mass} {cell.get_fluid(0).get_component(1).mass} {cell.get_fluid(0).get_component(2).mass} '
                        f'{cell.get_fluid(0).get_component(3).mass} {cell.get_fluid(0).get_component(4).mass} '
                        f'{cell.get_fluid(1).mass} '
                        f'{cell.get_fluid(2).mass} '
                        f'{cell.get_fluid(3).mass} '
                        f'{cell.get_fluid(4).get_component(0).mass} {cell.get_fluid(4).get_component(1).mass} '
                        f'{temp} {pres}\n')

def mass(i, j=None):
    """
    input 
    i = fluid
    j = component
    
    output 
    contour plot saturation
    """    
    X = []
    Z = [] 
    total = []
    fluid = []
    for cell in model.cells:
        x, y, z = cell.pos
        X.append(x)
        Z.append(z)
        total.append(cell.fluid_vol)
        if cell.get_fluid(i).component_number == 0 :
            fluid.append(cell.get_fluid(i).vol)
        else:
            fluid.append(cell.get_fluid(i).get_component(j).vol)
    
    
    vol_frac = [ i / j for i, j in zip(fluid, total)]
    vol_frac = vol_frac[:len(vol_frac) - 1]
    vol_frac = np.array(vol_frac)
    vol_frac = np.transpose(vol_frac.reshape(200, 200))
    fig, ax = pyplot.subplots()
    plot = ax.contourf(vol_frac, 20, extent=[0, 100, 0, 100], cmap='coolwarm', antialiased=True)
    ax.set_xlabel('x, m')
    ax.set_ylabel('y, m') 
    ax.set_ylim(100, 0)
    cbar = fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap), ax=ax)
    cbar.set_label('vol_sat')  # Set label for the colorbar
    
def temperature():
    """
    input 
    i = fluid
    j = component
    
    output 
    contour plot saturation
    """    
    X = []
    Z = [] 
    temp = []
    for cell in model.cells:
        x, y, z = cell.pos
        temp.append(cell.get_attr(seepage.cell_keys(model).temperature))
    
    
    temp = temp[:len(temp) - 1]
    temp = np.array(temp)
    temp = np.transpose(temp.reshape(200, 200))
    fig, ax = pyplot.subplots()
    plot = ax.contourf(temp, 20, extent=[0, 100, 0, 100], cmap='coolwarm', antialiased=True)
    ax.set_xlabel('x, m')
    ax.set_ylabel('y, m') 
    ax.set_ylim(100, 0)
    cbar = fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap), ax=ax)
    cbar.set_label('Temperature, K')  # Set label for the colorbar
    pyplot.show()

# temperature()

def iterate():
    name = os.path.basename(__file__)
    
    data_path = os.path.join(os.getcwd(), f'data_{name}')
    resu_path = os.path.join(data_path, f'Results_rate')
    os.makedirs(os.path.join(os.getcwd(), data_path, f'Results_rate'), exist_ok=True)
    folder = os.path.join(os.getcwd(), data_path, f'Results_rate')
    
    solver = ConjugateGradientSolver(tolerance=1.0e-8)
    for step in range(100000000):
        seepage.iterate(model, solver=solver)
        pre_ctrl.update(seepage.get_time(model))
        monitor.update(dt=seepage.get_dt(model))
        
        
        time = seepage.get_time(model) / (3600 * 24)
        if time > 3600 * 24 * 365 * 20:
            print(f'time Finish = {seepage.get_time(model) / 3600*34*365}')
            break
        
        if step % 1 == 0:
            temperature()
            monitor.save(os.path.join(folder, f'prod_{name}.txt'))
            print(f'time = {time}')
            
        
iterate()        



















    

