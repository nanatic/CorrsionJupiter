from typing import List

import json
from thermopack.cubic import cubic
from fluids.two_phase import Taitel_Dukler_regime

from PVT import *


from typing import List
import math
import pandas as pd
from thermopack.cubic import cubic
import copy
import sys
import time
import json

sys.path.insert(0, '../pycThermopack/')

from PVT import *
class Tube_point:
    name: str = 'mixture'  # name: string, name of the material
    phase_name: List[str] = ['gas']  # phase_name: list of strings, name of phases
    number_of_fluids: int = 1  # number_of_fluids: integer, number of phases
    temperature: float = 900.0  # T: double, local temperature, Kelvin
    pressure: float = 400325.0  # p: double, local pressure, Pascal
    velocity: float = 0.5  # v: double, mixture velocity, m/s
    diameter: float = 0.5  # D: double, tube local diameter, m
    length: float = 20  # length: double, tube length, m. Zero if it's endpoint
    section_type = 1  # section type: integer, type of the section for losses calculation, not used yet - enum!
    molar_composition: List[float] = [1.0]  # molar_composition: list of doubles, molar composition [probably constant]
    molar_masses: List[float] = [0.01604, 0.018015]  # molar_masses: list of doubles, molar masses [constant overall]
    vapor_components: List[float] = [
        0.5]  # vapor_components: list of doubles, vapor distribution over components (sum is 1)
    liquid_components: List[float] = [
        0.5]  # liquid_components: list of doubles, liquid distribution over components (sum is 1)
    components_density: List[float] = [ 665,997]  # components_density: list of doubles, density distribution over components
    overall_density: float = 100.0  # overall_density: double, overall density, kg/m^3
    overall_vapor_fraction: float = 0.1  # overall_vapor_fraction: double, vapor distribution over mixture
    overall_liquid_fraction: float = 0.9  # overall_liquid_fraction: double, liquid distribution over mixture
    liquid_viscosities: List[float] = [
        1.0e-2]  # liquid_viscosities: list of doubles, viscosity of liquid parts over components
    vapor_viscosities: List[float] = [
        1.0e-1]  # vapor_viscosities:  list of doubles, viscosity of vapor parts over components
    liquid_overall_viscosity: float = 1.0e-2  # liquid_overall_viscosity: double, viscosity of liquid part
    vapor_overall_viscosity: float = 1.0e-1  # vapor_overall_viscosity: double, viscosity of vapor part
    overall_viscosity: float = 1.0e-1  # overall_viscosity: double, viscosity of mixture
    flow_mode: str = "bubble"  # flow_mode: string, name of selected flow flow_mode
    flow_mode_key: float = 1.0  # flow_mode_key: double, currently XTT, later - other number to characterize flow_mode
    flow_mode_friction_factor: float = 1.0  # flow_mode_friction_factor: double, currently from XTT
    reynolds_number: float = 10000.0  # reynolds_number: double, Reynolds number for ...
    roughness = 0.01  # roughness: шероховатость внутренней поверхности трубы
    mass = 0.1
    angle = 0  # angle: угол наклона трубы
    list_length = []
    list_diameter = []
    components = ''


def start_point_from_excel(point):  # initialization of start point

    print(f'Ввод входных данных. в качестве десятичного разделителя используйте  "." (точку)')

    with open('Input_data.txt') as json_file:
        data = json.load(json_file)

    data['temperature'] = float(input("Введите температуру"))
    data['pressure'] = float(input("Введите давление"))

    list_ = []
    sub_molA = float(input(r"Введите молярную долю газа"))
    list_.append(sub_molA)
    sub_molB = float(input("Введите молярную долю жидкости"))
    list_.append(sub_molB)

    data['molar_composition'] = list_

    list_ = []
    sub_masA = float(
        input(r"Введите молярную массу газа"))
    list_.append(sub_masA)
    sub_masB = float(
        input("Введите молярную массу жидкости"))
    list_.append(sub_masB)

    data['molar_masses'] = list_
    data['diameter'] = float(input("Введите диаметр"))
    vg = float(input(f'Введите объемный расход газа за сутки'))
    vf = float(input(f'Введите объемный расход жидкости за сутки'))
    v = ((vg + vf) / 86400 / (math.pi * (data['diameter'])**2 / 4))

    data['velocity'] = v

    data['length'] = float(input("Введите длину трубопровода"))
    data['vapor_viscosities'] = float(input("Введите вязкость газа"))
    data['liquid_viscosities'] = float(input("Введите вязкость жидкости"))

    list_ = []
    sub_A = float(
        input(r"Введите плотность газа"))
    list_.append(sub_masA)
    sub_B = float(
        input("Введите плотность жидкости"))
    list_.append(sub_masB)
    data['components_density'] = list_

    data['roughness'] = float(input("Введите шероховатость стенок трубопровода [м] "))
    data['mass'] = float(input("Введите расход флюида (кг / с) "))
    data['angle'] = float(input("Введите угол наклона трубопровода"))



    lst = []
    dst = []
    n = int(input("Введите колличество отрезков трубопровода : "))

    for i in range(1, n+1):
        ele = int(input(f'Введите длину отрезка номер {i}'))
        diam = float(input(f'Введите диаметр отрезка номер {i}'))
        # adding the element
        lst.append(ele)
        dst.append(diam)
    data['list_length'] = lst
    data['list_diameter'] = dst

    # with open('sub_input.txt','w') as inpt:
    #     inpt.write(data)


    point.temperature = data['temperature']
    point.pressure = data['pressure']
    point.molar_composition = data['molar_composition']
    point.molar_masses = data['molar_masses']
    point.velocity = data['velocity']
    point.diameter = data['diameter']
    point.length = data['length']
    point.vapor_viscosities = data['vapor_viscosities']
    point.liquid_viscosities = data['liquid_viscosities']
    point.components_density = data['components_density']
    point.roughness = data['roughness']
    point.mass = data['mass']
    point.angle = data['angle']
    point.list_length = data['list_length']
    point.list_diameter = data['list_diameter']
    point.components = data['components']

    update_point_state(point)


def define_tube_params(point, diameter, length, density_old):
    q = point.velocity * point.diameter * point.diameter * density_old
    new_velocity = q / (diameter * diameter * point.overall_density)  # mass balance, pi/4 is skipped
    # due to presence in both parts of equation
    point.diameter = diameter
    point.length = length
    point.velocity = new_velocity


def update_point_state(point):
    """
    Updates tube point parameters after changing local temperature and pressure
    """
    rk_fluid = cubic(point.components, 'SRK')  # obsolete
    x, y, vap_frac, liq_frac, phase_key = rk_fluid.two_phase_tpflash(point.temperature, point.pressure,
                                                                     point.molar_composition)
    point.vapor_components = x
    point.liquid_components = y
    point.overall_vapor_fraction = vap_frac
    point.overall_liquid_fraction = liq_frac

    temp, = rk_fluid.specific_volume(point.temperature, point.pressure, point.molar_composition, 1)
    density_1 = point.molar_masses[0] / temp

    temp, = rk_fluid.specific_volume(point.temperature, point.pressure, point.molar_composition, 2)
    density_2 = point.molar_masses[1] / temp

    point.components_density = [density_1, density_2]
    point.overall_density = calculate_overall_density(point)
    ethanol_viscosity = ethanol_viscosity_from_temperature(point.temperature)
    n2_viscosity = n2_viscosity_from_temperature(point.temperature)
    point.liquid_viscosities = [ethanol_viscosity, n2_viscosity]
    point.vapor_viscosities = [ethanol_viscosity, n2_viscosity]


def calculate_Re(point):
    """
    Calculates the Reynolds number based on the total density and total viscosity of the medium.
    :return: Reynolds number
    """
    return point.velocity * point.diameter * point.overall_density / point.overall_viscosity


def calculate_xtt(point):
    """
    Calculates the parameter by which the flow mode can be obtained.
    NOTE: The simplest correlation has been applied, which will require adjustments in the future
    :return: xtt - Lockhart-Martinelli parameter
    """
    liquid_density = point.components_density[0]
    gas_density = point.components_density[1]
    liquid_viscosity = point.liquid_viscosities[0]  # ? liquid_overall_viscosity?
    gas_viscosity = point.liquid_viscosities[1]  # ?
    velocity = point.velocity
    diameter = point.diameter
    return ((1.096 / liquid_density) ** 0.5) * ((liquid_density / gas_density) ** 0.25) * (
            (gas_viscosity / liquid_viscosity) ** 0.1) * ((velocity / diameter) ** 0.5)


def calculate_overall_density(point):  # необходимо дописать учёт агрегатного состояния
    return sum(point.molar_composition[i] * point.components_density[i] for i in range(point.number_of_fluids))


def calculate_lambda(point):
    if point.reynolds_number < 2300:
        return 64 / point.reynolds_number
    else:
        return 0.316 / (point.reynolds_number ** 0.25)


def calculate_pressure_loss(point):
    xi = calculate_lambda(point) * point.length / point.diameter
    return (xi * point.velocity ** 2) * 0.5 * point.overall_density


def calculate_viscosity(point, friction_factor):
    liquid_viscosity = point.liquid_viscosities[0]  # ? liquid_overall_viscosity?
    gas_viscosity = point.liquid_viscosities[1]  # ?
    return friction_factor * liquid_viscosity + (1 - friction_factor) * gas_viscosity



def return_mode(xtt):
    list_ = []
    list_ = Taitel_Dukler_regime(m=Tube_point.mass, x=0.7, rhol=Tube_point.components_density[1],
                                 rhog=Tube_point.components_density[0], mul=Tube_point.liquid_viscosities[0],
                                 mug=Tube_point.vapor_viscosities[0], D=Tube_point.diameter, angle=Tube_point.angle,
                                 roughness=Tube_point.roughness)
    return list_[0]


# def return_mode(xtt):
#     if xtt < 10: return 'bubble'
#     if 10 <= xtt < 100:
#         return 'plug'
#     if 100 <= xtt < 1000:
#         return 'slug'
#     if 1000 <= xtt < 10000:
#         return 'annular'
#     if 10000 <= xtt:
#         return 'mist'
#     return 'undefined'


# liquid to solid viscosity calculation:

def return_friction_factor(xtt):
    """
    Outputs the friction factor to calculate the viscosity.
    :param xtt:
    :return:
    """
    if xtt < 10:
        return 1
    if 10 <= xtt < 100:
        return 0.9
    if 100 <= xtt < 1000:
        return 0.8
    if 1000 <= xtt < 10000:
        return 0.7
    if 10000 <= xtt:
        return 0.6
    return 0

def pvt_block(point: Tube_point, new_pressure, new_temperature):
    point.temperature = new_temperature
    point.pressure = new_pressure
    update_point_state(point)
    return point
