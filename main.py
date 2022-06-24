import numpy as np
import cantera as ct
import math
import csv
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.yaml')

#warunki poczatkowe - podstawowe zmienne - temperatura i ciśnienie
init_T = 350.0 #K
init_P = 1 #atm

#parametry w zbiorniku paliwa - propan
fuel_T = init_T #temperatura
fuel_P = init_P*ct.one_atm*50 #cisnienie wyrażone w Pascalach
fuel_X = 'C3H8:1.0' #ilosc moli
gas.TPX = fuel_T, fuel_P, fuel_X
fuel = ct.Reservoir(gas)
fuel_mw = gas.mean_molecular_weight
fuel_k = gas.cp/gas.cv


#parametry w zbiorniku utleniacza
oxidizer_T = init_T
oxidizer_P = init_P*ct.one_atm*50
oxidizer_X = 'O2:5.0'
gas.TPX = oxidizer_T, oxidizer_P, oxidizer_X
oxidizer = ct.Reservoir(gas)
oxidizer_mw = gas.mean_molecular_weight
oxidizer_k = gas.cp/gas.cv

#zrodlo zaplonu - wolne rodniki wodoru
gas.TPX = 300.0, ct.one_atm, 'H:1.0'
igniter = ct.Reservoir(gas)

#komora spalania
gas.TPX = 300.0, 1.1*ct.one_atm, 'O2:5.0'   #komora spalania jest poczatkowo wypelniona tlenem
cchamber = ct.IdealGasReactor(gas)
cchamber.volume = 0.015 #m^3
cchamberP = gas.P

#dysza
gas.TPX = 300.0, 1*ct.one_atm, 'O2:5.0'
exhaust = ct.Reservoir(gas)
exaustP = ct.one_atm

#zdefiniowanie potrzebnych funkcji
def kappa(gas):
    return gas.cp/gas.cv

def critical_flow(gasin, gasinP, gasinT, gasinmw, k, gasoutP, area):
    R = ct.gas_constant/gasinmw
    return (area*gasinP*math.sqrt(k/(R*gasinT))*(2/(k+1))**((k+1)/(2*(k-1))))/(gasinP - gasoutP)

v1 = ct.Valve(fuel, cchamber)
v2 = ct.Valve(oxidizer, cchamber)
v3 = ct.Valve(cchamber, exhaust)

#Konfiguracja zapalnika
fwhm = 0.008
amplitude = 0.01
t0 = 0.05
igniter_mdot = lambda t: amplitude * math.exp(-(t - t0) ** 2 * 4 * math.log(2) / fwhm ** 2)
m3 = ct.MassFlowController(igniter, cchamber, mdot=igniter_mdot)


states = ct.SolutionArray(gas, extra=['t','vel','thrust'])

#Symulacja zawiera tylko jeden reaktor - komora splania
sim = ct.ReactorNet([cchamber])
time = 0.0
tfinal = 0.01



#Petla obliczeniowa
while time < tfinal:

    v1.valve_coeff = 4e-5
    v2.valve_coeff = 4e-5
    v3.valve_coeff = 5e-4
    time = sim.step()

    v = (2 * (kappa(gas) * ct.gas_constant / gas.mean_molecular_weight) / (kappa(gas) - 1)
    * gas.T * (1 - (exaustP / gas.P)) ** ((kappa(gas) - 1) / kappa(gas))) ** 0.5
    T = gas.density * v3.valve_coeff * v * v
    states.append(gas.state, t=time, vel=v, thrust=T)
    print(round(time, 4), round(gas.P / 1e6, 3), round(gas.T,2), round(gas.density,1), round(v,2))



#Plotowanie wykresu predkosci wylotowej, zmieniajacej sie w czasie
plt.xlabel('$t[s]$')
plt.ylabel('$Thrust [N]$')
plt.plot(states.t, states.thrust)
plt.show()

#Plotowanie wykresu ciągu
plt.xlabel('$t[s]$')
plt.ylabel('$V [m/s]$')
plt.plot(states.t, states.vel)
plt.show()


