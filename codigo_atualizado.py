import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect

def loitter(E, C, LD):
    return np.exp(-E*C/LD)

def takeoff():
    return 0.97

def landing():
    return  0.995

def climb():
    return 0.985

def climb_acelerado(M):
    return 0.991 - 0.00325*M

def cruise(R, C, LD, V):
    return np.exp(-R*C/V/LD)


def alijamento(R, C, LD, V, m_seed):
    return lambda Wi: Wi*np.exp(-R*C/(V*LD)) - (1 - np.exp(-R*C/(V*LD)))*(LD/C)*m_seed*g


if __name__=='__main__':
    LDmax =     12
    LD_loitter = 0.886*LDmax
    E_loitter = 20*60        # s
    C_cruise =  0.5/3600     # 1/s
    C_loitter = 0.4/3600     # 1/s
    mseed =     5/1000       #vazão mássica de sementes em kg/s
    mass_seed = 63           # kg
    t_aerial_seeding = mass_seed/mseed # s
    V_aerial_seeding = 20    # m/s
    V_cruise = 30            # m/s (Mach ~ 0.1) 
    R_aerial_seeding = V_aerial_seeding*t_aerial_seeding # m
    print(f'Distância de alijamento {R_aerial_seeding/1000} [km]')
    print(f'Tempo de alijamento {t_aerial_seeding/3600} [h]')
    R_cruise = 1000          # m
    g =         9.81         # m/s^2
    Wpl =       70*g         # N
    
    massa_combustível = 12 # kg

    W1_W0 = takeoff()

    W2_W1 = climb_acelerado(M=0.1)

    W3_W2 = cruise(R_cruise, C_cruise, LDmax, V_cruise)

    W4 = alijamento(R_aerial_seeding, C_cruise, LDmax, V_aerial_seeding, mseed)

    W5_W4 = climb_acelerado(M=0.1)

    W6_W5 = cruise(R_cruise, C_cruise, LDmax, V_cruise)

    W7_W6 = landing()

    
    # Na equação abaixo We = 0.49*W0**0.96  
    f_W0 = lambda W0 : 0.49*W0**0.96 + 1.06*(W0 - W7_W6*W6_W5*W5_W4*W4(W0*W1_W0*W2_W1*W3_W2) - mass_seed*g) - W0 + Wpl


    # x = np.linspace(0, 150*g, 10000)
    # plt.plot(x, f_W0(x))
    # plt.xlabel(r'$W_0$ [N]')
    # plt.ylabel('Resíduo')
    # plt.show()

    W0 = bisect(f_W0, 0, 150*g)


    print(f'MTOW {W0} [N]')
    print(f'MTOM {W0/9.81} [kg]')
    print(f'W_fuel {1.06*(W0 - W7_W6*W6_W5*W5_W4*W4(W0*W1_W0*W2_W1)) - (1.06*mass_seed*g)} [N]')
    print(f'M_fuel {(1.06*(W0 - W7_W6*W6_W5*W5_W4*W4(W0*W1_W0*W2_W1)) - (1.06*mass_seed*g))/9.81} [kg]')
    print(f'Peso vazio da aeronave {0.49*W0**0.9} [N]')
    print(f'Carga paga {Wpl/g} [kg]')
    print(f'Massa de sementes {mass_seed} [kg] \n' )


    print(f'Takeoff: {W1_W0}')
    print(f'Climb: {W2_W1}')
    print(f'Cruise 1: {W3_W2}')
    print(f'Aerial Seeding: {W4(W0*W1_W0*W2_W1*W3_W2)/(W0*W1_W0*W2_W1)}')
    print(f'Climb 2: {W5_W4}')
    print(f'Cruise 2: {W6_W5}')
    print(f'Landing: {W7_W6}')


   


