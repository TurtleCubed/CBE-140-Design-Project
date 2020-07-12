import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

RUN_TYPE = "-"

# I made three modes of running this code so I could test various parts of this problem
# "single" tests one set of parameters. I used this to test if my math was correct
if RUN_TYPE is "single":
    FLOW_IN = 100  # kg per second
    REACTOR_DIAMETER = 10  # meters
    REACTOR_TEMPERATURE = 300  # Celsius
# "graph" tests various diameters and temperatures. This allowed me to visualize where the maximum would be. I found
# that the maximum would fall on one of the edges of the range since there was no relative maxima.
elif RUN_TYPE is "graph":
    FLOW_IN = 399  # kg per second
    REACTOR_DIAMETER = np.linspace(4, 10, 20)  # meters
    REACTOR_TEMPERATURE = np.linspace(25, 600, 500)  # Celsius
# This case is what I used to optimize profit. It varies flow rate, reactor size, and reactor temperature.
else:
    FLOW_IN = np.linspace(10, 1000, 10000)  # kg per second
    REACTOR_DIAMETER = np.linspace(4, 10, 5)  # meters
    REACTOR_TEMPERATURE = np.linspace(25, 600, 10)  # Celsius

# Since our project statement assumes the reaction follows some equilibrium, I set a minimum amount of reaction time
# to allow the reaction to proceed to equilibrium.
MINIMUM_RESIDENCE_TIME = 3600  # seconds

# These are various specifications set forth in the prompt that are useful in calculating profit.
m_cost = 0.04 * 1000  # $ per kilogram
p_cost = 0.12 * 1000  # $ per kilogram
p_molmass = 100 * 0.001  # kilograms per mole
m_molmass = p_molmass / 2  # kilograms per mole
m_cp = 100 * 1000  # joules per kilogram per Kelvin
m_density = 822  # kg per m^3
p_density = 1030  # kg per m^3 ASSUMPTION: M dimerizes like ethylene oxide into dioxane
h_rxn = 20 * 1000  # joules per mole
h_t_c = 2.5  # watts per m^2 per Kelvin
r_gas = 8.314  # joules per mol per kelvin

# ASSUMPTION: the heat capacity for M and P are the same
if True:  # set to False for linear heat capacity scaling
    p_cp = m_cp
else:
    p_cp = m_cp / m_molmass * p_molmass


# Uses the given equilibrium to solve for the output flow rate of P, see write-up for derivation
def concentration_p(conc_m_in, k):
    x = -1 * (np.sqrt(8 * k * conc_m_in + 1) - 4 * k * conc_m_in - 1) / (8 * k)
    return x


# Calculates for the mass ratio of P to M in the output stream
def output_ratio(conc_m_in, conc_p):
    num = conc_p * p_molmass
    den = (conc_m_in - (2 * conc_p)) * m_molmass
    return num / den


# Calculates the output mass flowrate based on the output molar ratio
def output_composition(flow_rate, ratio):
    mp = flow_rate * ratio / (ratio + 1)
    mm = flow_rate / (ratio + 1)
    return mp, mm


# Uses the given k statement to calculate for k
def k_rxn(reactor_temperature):
    return np.exp(-1 * h_rxn / (r_gas * (273 + reactor_temperature)))


# Uses geometry to calculate for reactor volume
def reactor_volume(diameter):
    return np.pi * (diameter / 2) * (diameter / 2) * 2 * diameter


# Uses geometry to calculate for reactor area. Used in determining heat loss and replacement cost
def reactor_area(diameter):
    height = 2 * diameter
    base_area = np.pi * (diameter / 2) * (diameter / 2)
    wall_area = np.pi * diameter * height
    area = base_area + wall_area
    return area


# Uses the given material cost to calculate for replacement costs
def reactor_cost(diameter):
    material_cost = 150000  # dollars per m^2
    cost_per_year = material_cost * reactor_area(diameter)
    return cost_per_year / 365  # dollars per day


# Calculates for the required heat input to the reactor. Net heat is calculated from a few enthalpies: the heat lost
# from the surface of the reactor, the heat absorbed by the reaction, and the heat required to raise M and P from
# the standard temperature to the specified temperature.
def heat_in(flow_in, p_out, zeta, reactor_temp erature, diameter):
    loss = 2.5 * reactor_area(diameter) * (reactor_temperature - 25)
    delta_h_rxn = h_rxn * zeta
    delta_h_m = flow_in * m_cp * (reactor_temperature - 25)
    delta_h_p = p_out * p_cp * (reactor_temperature - 25)
    return loss + delta_h_rxn + delta_h_m + delta_h_p  # joules per second


# Calculates for profit. Money is earned through selling P. Money is lost through purchasing M, heating the reactor,
# and replacing the reactor material.
def profit(flow_in, diameter, reactor_temperature):
    conc_m_in = m_density / (1000 * m_molmass)  # Pure M has an intrinsic concentration
    k = k_rxn(reactor_temperature)
    # Using the above functions to calculate for output flow rates
    mass_p_out, mass_m_out = output_composition(flow_in, output_ratio(conc_m_in, concentration_p(conc_m_in, k)))
    mol_p_out = mass_p_out / p_molmass
    # Buying cost is a function of the flow rate
    m_buying_cost = m_cost * flow_in * 3600 * 24  # dollars per day
    p_revenue = mass_p_out * p_cost * 3600 * 24  # dollars per day
    rxr_cost = reactor_cost(diameter)  # dollars per day
    daily_energy_usage = heat_in(flow_in, mass_p_out, mol_p_out,
                                 reactor_temperature, diameter) * 3600 * 24  # joules per day
    daily_energy_usage = daily_energy_usage / 0.12  # 12% heat input efficiency
    energy_cost = daily_energy_usage / np.power(10, 9) * 14  # dollars per day
    # Calculates residence time using flow rates
    tau = reactor_volume(diameter) / (mass_p_out / p_density + mass_m_out / m_density)
    return p_revenue - rxr_cost - energy_cost - m_buying_cost, tau


# Used for testing functionality
if RUN_TYPE is "single":
    print(profit(FLOW_IN, REACTOR_DIAMETER, REACTOR_TEMPERATURE))

# Used for testing ranges
elif RUN_TYPE is "graph":
    fig = plt.figure(figsize=plt.figaspect(0.5))

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.set_xlabel('Reactor Diameter (m)')
    ax.set_ylabel('Reactor Temperature (C)')
    ax.set_zlabel('Profit ($ / day)')

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.set_xlabel('Reactor Diameter (m)')
    ax2.set_ylabel('Reactor Temperature (C)')
    ax2.set_zlabel('Residence Time (s)')

    X, Y = np.meshgrid(REACTOR_DIAMETER, REACTOR_TEMPERATURE)
    Z, tau = profit(FLOW_IN, X, Y)

    for i in range(len(REACTOR_DIAMETER)):
        for j in range(len(REACTOR_TEMPERATURE)):
            if tau[j][i] < MINIMUM_RESIDENCE_TIME:
                Z[j][i] = 0

    ax.plot_surface(X, Y, Z)
    ax2.plot_surface(X, Y, tau)
    plt.show()

# Used to optimize profit for a wide range of values of all three independent varibles
else:
    big_table = np.zeros((len(FLOW_IN), len(REACTOR_DIAMETER), len(REACTOR_TEMPERATURE)))
    max_profit = 0
    optimal_config = [0, 0, 0]
    for i in range(len(FLOW_IN)):
        for j in range(len(REACTOR_DIAMETER)):
            for k in range(len(REACTOR_TEMPERATURE)):
                big_table[i, j, k], tau = profit(FLOW_IN[i], REACTOR_DIAMETER[j], REACTOR_TEMPERATURE[k])
                if tau < MINIMUM_RESIDENCE_TIME:
                    big_table[i, j, k] = 0
                if big_table[i, j, k] > max_profit:
                    max_profit = big_table[i, j, k]
                    optimal_config = [i, j, k]
    best_flow = FLOW_IN[optimal_config[0]]
    best_diam = REACTOR_DIAMETER[optimal_config[1]]
    best_temp = REACTOR_TEMPERATURE[optimal_config[2]]
    print("Maximum profit of $", max_profit, " per day achieved with:")
    print(best_flow, " kg M per second")
    print(best_diam, " m diameter tank")
    print(best_temp, " °C")
