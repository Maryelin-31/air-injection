# -*- coding: utf-8 -*-
from zml import *
from zmlx.alg import make_fname

import matplotlib.tri as tri
from matplotlib.cm import ScalarMappable
from matplotlib import pyplot as pyplot
import numpy as np
import os

from mesh import *
from Reaction_conf import create
from initial_cond import *

from zmlx.kr.create_krf import create_krf
from zmlx.alg.linspace import linspace
from zmlx.utility.HeatInjector import HeatInjector
from zmlx.utility.RelativePermeability import stone_model_II

"following Profesor Workspace"

def seepage_script(rate_inj):

    "Create Mesh"
    mesh = create_mesh()
    # print(mesh)
    
    config = create()
    
    "Initial Conditions Reservoir"
    ini = create_initial()
    
    "create model"
    model = config.create(mesh, **ini)  # create
    model.set(gravity=(0, -10, 0))
    
    "relative permeability"
    model.set_kr(config.iGAS, *create_krf(0.05, 3.0))
    model.set_kr(config.iwat, *create_krf(0.05, 3.0))
    model.set_kr(config.iLO, *create_krf(0.05, 3.0))
    model.set_kr(config.iHO, *create_krf(0.05, 3.0))
    
    "Production Well"
    
    pos = (15.0, 8.0, 1000)
    cell = model.get_nearest_cell(pos=pos)
    virtual_cell = config.add_cell(model, vol=cell.get_attr(config.cell_keys['vol']), pos=pos, porosity=1.0e3, pore_modulus=20e6,
                                   temperature=ini['temperature'](*pos), p=1.0e6,
                                   s=ini['s'](*pos))
    x, y, z = cell.pos
    virtual_cell.set_attr(config.cell_keys['g_heat'], 1.0e-20)
    face = config.add_face(model, virtual_cell, model.get_cell(cell.index),
                           heat_cond=0, perm=1e-14, area=1.0, length=1.0, )
   
    
    prectrl = PressureController(virtual_cell, t=[0, 1e10], p=[1e6, 1e6])
    monitor = SeepageCellMonitor(get_t=lambda: config.get_time(model),
                                  cell=(virtual_cell, prectrl))
    
    "Injection Configuration"
    cell_ids = set()
    z = 0
    y = 5.0
    for x in (10.0, 20.0):
        # for y in (5.0, 10.0):
        cell = model.get_nearest_cell(pos=(x, y, z))
        # print(cell)
        cell_ids.add(cell.index)
    
       
    # Add_injector
    for cell_id in cell_ids:
        fa_t = 0
        fa_c = 1
        cell1 = model.get_cell(cell_id)
        flu = cell1.get_fluid(0).get_component(3)  # Air
        flu.set_attr(fa_t, 700)
        flu.set_attr(fa_c, 1000)
        flu.set_attr(config.flu_keys['temperature'], 700)
        model.add_injector(cell=cell1, fluid_id=[0, 3], flu=flu,
                            pos=cell1.pos, radi=1.0, opers=[(0, rate_inj)])
        

    
    "Time step Strategy"
    config.set_dv_relative(model, 0.5)  # The ratio of the distance traveled by each time step to the grid size
    config.set_dt(model, 0.01)  # initial value for time step
    config.set_dt_max(model, 24 * 3600 * 7)  # Maximum value of time step <one week> # Maximum value of time step <one week>
    config.set_dt_min(model, 3600)  # Minimum step size is 1 hour
    
    "Plots"
    
    def Temperature_cell(time, step):  # Matrix
    
        custom_colors = [
            [0, 0, 1], [0, 0.1875, 1], [0, 0.375, 1], [0, 0.5625, 1],
            [0, 0.75, 1], [0, 0.9375, 1], [0.125, 1, 0.875], [0.3125, 1, 0.6875],
            [0.5, 1, 0.5], [0.6875, 1, 0.3125], [0.875, 1, 0.125], [1, 0.9375, 0],
            [1, 0.75, 0], [1, 0.5625, 0], [1, 0.375, 0], [1, 0.1875, 0]]
        custom_cmap = pyplot.cm.colors.ListedColormap(custom_colors)
    
        name = os.path.basename(__file__)
        X = []
        Y = []
        Temp = []
        i = 120
        j = 60
        for cell in model.cells:
            x, y, z = cell.pos
            X.append(x)
            Y.append(y)
            Temp.append(cell.get_attr(config.cell_keys['temperature']))
    
        if len(Temp) == mesh.cell_number:
            Temp = Temp
        else:
            Temp = Temp[:len(Temp) - 1]
    
        # Temp = Temp[:len(Temp) - 1]
        [xx, yy] = np.meshgrid(X, Y)
        fig, ax = pyplot.subplots()
        Temp = np.array(Temp)
        Temp = np.transpose(Temp.reshape(i, j))
        tini = 300
        tmax = 1000  # temperature
        plot = pyplot.contourf(Temp, 20, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=custom_cmap, vmin=np.min(tini),
                               vmax=np.max(tmax))
        Time = round(time / (3600 * 24), 1)
        pyplot.title(f'Temperature of Rock (K) at {Time} days', fontsize=10)
        pyplot.xlabel('x/m', fontsize=10)
        pyplot.ylabel('y/m', fontsize=10)
        pyplot.tick_params(axis='both', which='major', labelsize=10)
        pyplot.ylim(max(Y), min(Y))
        pyplot.clim(vmin=np.min(tini), vmax=np.max(tmax))
        fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap))
    
    
        # PlotPath = f'Temperature_{step}.png'   
        # folder = os.path.join(os.getcwd(),f'Temp_cell_Plot_{name}')
        # plotfile = os.path.join(folder, PlotPath)
        # plt.savefig(plotfile , format='png')
        pyplot.show()
        # plt.close()
    
    # Temperature_cell(0, 0)
    
    def Temperature_fluid(time, step, fid, i):  # Matrix
    
        custom_colors = [
            [0, 0, 1], [0, 0.1875, 1], [0, 0.375, 1], [0, 0.5625, 1],
            [0, 0.75, 1], [0, 0.9375, 1], [0.125, 1, 0.875], [0.3125, 1, 0.6875],
            [0.5, 1, 0.5], [0.6875, 1, 0.3125], [0.875, 1, 0.125], [1, 0.9375, 0],
            [1, 0.75, 0], [1, 0.5625, 0], [1, 0.375, 0], [1, 0.1875, 0]]
        custom_cmap = pyplot.cm.colors.ListedColormap(custom_colors)
        
        name = os.path.basename(__file__)
        X = []
        Y = []
        Temp = []
        i = 120
        j = 60
        for cell in model.cells:
            x, y, z = cell.pos
            X.append(x)
            Y.append(y)
            if cell.get_fluid(fid).component_number == 0:
                flu = cell.get_fluid(fid)
                Temp.append(flu.get_attr(config.flu_keys['temperature']))
            else:
                flu = cell.get_fluid(fid).get_component(i)
                Temp.append(flu.get_attr(config.flu_keys['temperature']))
    
        if len(Temp) == mesh.cell_number:
            Temp = Temp
        else:
            Temp = Temp[:len(Temp) - 1]
    
        [xx, yy] = np.meshgrid(X, Y)
        fig, ax = pyplot.subplots()
        Temp = np.array(Temp)
        Temp = np.transpose(Temp.reshape(i, j))
        tini = 300
        tmax = 1000  # temperature
        plot = pyplot.contourf(Temp, 20, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=custom_cmap, vmin=np.min(tini),
                            vmax=np.max(tmax))
        Time = round(time / (3600 * 24), 1)
        pyplot.title(f'Temperature of Fluid (K) at {Time} days', fontsize=10)
        pyplot.xlabel('x/m')
        pyplot.ylabel('y/m')
        pyplot.ylim(max(Y), min(Y))
        pyplot.clim(vmin=np.min(tini), vmax=np.max(tmax))
        fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap))
    
        # PlotPath = f'Temperature_{step}.png'   
        # folder = os.path.join(os.getcwd(),f'Temp_cell_Plot_{name}')
        # plotfile = os.path.join(folder, PlotPath)
        # plt.savefig(plotfile , format='png')
        pyplot.show()
        # plt.close()
         
    # Temperature_fluid(0, 0, fid=4, i=0) 
    
    def Pressure_cell(time, step):  # Matrix
    
        custom_colors = [
            [0, 0, 1], [0, 0.1875, 1], [0, 0.375, 1], [0, 0.5625, 1],
            [0, 0.75, 1], [0, 0.9375, 1], [0.125, 1, 0.875], [0.3125, 1, 0.6875],
            [0.5, 1, 0.5], [0.6875, 1, 0.3125], [0.875, 1, 0.125], [1, 0.9375, 0],
            [1, 0.75, 0], [1, 0.5625, 0], [1, 0.375, 0], [1, 0.1875, 0]]
        custom_cmap = pyplot.cm.colors.ListedColormap(custom_colors)
    
        name = os.path.basename(__file__)
        X = []
        Y = []
        press = []
        i = 120
        j = 60
        for cell in model.cells:
            x, y, z = cell.pos
            X.append(x)
            Y.append(y)
            press.append(cell.get_attr(config.cell_keys['pre']))
    
        if len(press) == mesh.cell_number:
            press = press
        else:
            press = press[:len(press) - 1]
    
        # press = press[:len(press) - 1]
        [xx, yy] = np.meshgrid(X, Y)
        fig, ax = pyplot.subplots()
        press = np.array(press)
        press = np.transpose(press.reshape(i, j))
        # pmin = 3.5e6
        # pmax = 1.0e8
        plot = pyplot.contourf(press, 20, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=custom_cmap, vmin=np.min(press),
                            vmax=np.max(press))
        Time = round(time / (3600 * 24), 1)
        pyplot.title(f'Pressure of Rock (Pa) at {Time} days', fontsize=10)
        pyplot.xlabel('x/m')
        pyplot.ylabel('y/m')
        pyplot.ylim(max(Y), min(Y))
        pyplot.clim(vmin=np.min(press), vmax=np.max(press))
        fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap))
    
        # PlotPath = f'Temperature_{step}.png'   
        # folder = os.path.join(os.getcwd(),f'Temp_cell_Plot_{name}')
        # plotfile = os.path.join(folder, PlotPath)
        # plt.savefig(plotfile , format='png')
        pyplot.show()
        # plt.close()   
    
    def mass(time, fid, i=None):  # fid = fluid , i = component
    
        custom_colors = [
             [0, 0, 1], [0, 0.1875, 1], [0, 0.375, 1], [0, 0.5625, 1],
             [0, 0.75, 1], [0, 0.9375, 1], [0.125, 1, 0.875], [0.3125, 1, 0.6875],
             [0.5, 1, 0.5], [0.6875, 1, 0.3125], [0.875, 1, 0.125], [1, 0.9375, 0],
             [1, 0.75, 0], [1, 0.5625, 0], [1, 0.375, 0], [1, 0.1875, 0]]
        custom_cmap = pyplot.cm.colors.ListedColormap(custom_colors)
         
        mass = []
        X = []
        Y = []
        i = 120
        j = 60
        for cell in model.cells:
            x, y, z = cell.pos
            X.append(x)
            Y.append(y)
            if cell.get_fluid(fid).component_number == 0:
                mass.append(cell.get_fluid(fid).mass)
            else:
                mass.append(cell.get_fluid(fid).get_component(i).mass)
    
        if len(mass) == mesh.cell_number:
            mass = mass
        else:
            mass = mass[:len(mass) - 1]
    
        [xx, yy] = np.meshgrid(X, Y)
        fig, ax = pyplot.subplots()
        mass = np.array(mass)
        mass = np.transpose(mass.reshape(i, j))
        plot = pyplot.contourf(mass, 20, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=custom_cmap, vmin=np.min(mass),
                               vmax=np.max(mass))
        Time = round(time / (3600 * 24), 1)
        pyplot.title(f'Mass (Kg) at {Time} days', fontsize=10)
        pyplot.xlabel('x/m', fontsize=8)
        pyplot.ylabel('y/m', fontsize=8)
        pyplot.tick_params(axis='both', which='major', labelsize=8)
        pyplot.ylim(max(Y), min(Y))
        pyplot.clim(vmin=np.min(mass), vmax=np.max(mass))
        fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap))
        # PlotPath = f'Temperature_{step}.png'   
        # folder = os.path.join(os.getcwd(),f'Temp_cell_Plot_{name}')
        # plotfile = os.path.join(folder, PlotPath)
        # plt.savefig(plotfile , format='png')
        pyplot.show()
        # plt.close()
    
    
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
        result_folder = os.path.join(os.getcwd(), 'data', f'Results_rate_{rate_inj}', 'Saturation_cells')
        
        if os.path.exists(f'result_folder'):
            import shutil
            shutil.rmtree(f'result_folder')            
        os.makedirs(result_folder, exist_ok=True)
        
        SavePath = os.path.join(result_folder, path)
        with open(SavePath, 'w') as file:
            for cell in model.cells:
                x, y, z = cell.pos
                satu = fluid_sat(cell)
                satu_str = ' '.join(str(i) for i in satu)
                file.write(f'{x} {y} {z} {satu_str}\n')

    
    def save_mass(path):
        result_folder = os.path.join(os.getcwd(), 'data', f'Results_rate_{rate_inj}', 'Results_cells')
        
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
    
    
    def run():
        data_path = os.path.join(os.getcwd(), 'data')
        resu_path = os.path.join(data_path, f'Results_rate_{rate_inj}')
        os.makedirs(os.path.join(os.getcwd(), 'data', f'Results_rate_{rate_inj}'), exist_ok=True)

        # folder = os.getcwd()
        folder = os.path.join(os.getcwd(), 'data', f'Results_rate_{rate_inj}')
        solver = ConjugateGradientSolver()
        solver.set_tolerance(1.0e-13)
        
        
        " initial pression and permeability"
        prei = 10e-6
        kini = 1.0e-15
        cfr= 0.05e-6 #average of compresibility fracture

        for step in range(1000000):
            config.iterate(model, solver=solver)
            time = config.get_time(model)
            
            # if time > 3600 * 24 * 100:
            #     face_connect_well.cond = face_backup
            
            "Permeability vs Pressure"

            # for face in model.faces:
            #     pre = (face.get_cell(0).pre + face.get_cell(1).pre) / 2  
            #     dp = prei - pre
            #     perm= kini * (1 - cfr * (dp))**3 # Saidi (1987)
            #     face.set_attr(config.face_keys['perm'], perm)
            #     cond= perm * (face.get_attr(config.face_keys['area']) / face.get_attr(config.face_keys['length']))
            #     face.set_attr(config.face_keys['g0'], cond)
    
            # HeatCapacity Abdelagatov, 2021 base on Wittington, 2009
            for cell in model.cells:
                vol = cell.get_attr(config.cell_keys['vol'])
                temp = cell.get_attr(config.cell_keys['temperature'])
                heat_cap = (0.849315 + 0.000565 * temp - 25002.318124 * temp ** (-2)) * 1000
                cell.set_attr(config.cell_keys['mc'], vol * 2600 * heat_cap)
    
            # Heat Conductivity Jin et al, 2022 (experiment)
            ca_t = config.cell_keys['temperature']
            fa_g = config.face_keys['g_heat']
            for face in model.faces:
                temp = (face.get_cell(0).get_attr(ca_t) + face.get_cell(1).get_attr(ca_t)) / 2
                heat_con = (-0.0006 * temp) + 0.8126
                face.set_attr(fa_g, heat_con)
    
            if config.get_time(model) > 3600 * 24 * 365 * 3:
                print(f'{step}, Finish')
                break
            time = config.get_time(model)
            
    
            path = f'time_{time}.txt'
            if step % 100 == 0:
                save_mass(path)
                save_sat(path)
                monitor.update(dt=config.get_dt(model))
                # Temperature_cell(time, step)
                # Pressure_cell(time, step)
                # mass(time, fid=2, i=None)
                # mass(time, fid=3, i=None)
                monitor.save(os.path.join(folder, 'prod.txt'))
                print(f'{step} {time / (3600 * 24)}')
                    
    
    run()
        
    
seepage_script(rate_inj=72.5)
