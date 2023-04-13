#Dipole Matrix Element Calculator
#Written By: Frank Corapi

#This program calculates the dipole matrix element for a atomic transition
#which can then be used to calculate the expected Rabi frequencies. 

from sympy.physics.wigner import wigner_3j, wigner_6j
import numpy as np

#Constants
h = 6.626e-34 #Js (Planck's constant)
hbar = h/(2*np.pi) #Js (reduced Planck's constant)
e = 1.602e-19 #C (charge of electron)
a_0 = 5.29177210903e-11 #m (Bohr Radius)
e_0 = 8.85418782e-12 #(permittivity of free space in SI units)
c = 2.998e8 #m/s (speed of light)

#******Input Parameters************************************

#Light properties
k = 1 #Rank of the tensor operator (1 for dipole, 2 for quadrupole, etc..)
q = -1 #Polarization of light  (0 for pi, 1 for sigma+, -1 for sigma-)
gamma = 3.734e7 #2*np.pi*6.035e6 #s^-1 (Natural linewidth/Einstein A coefficient)
wavelength = 769.896456e-9 #770.108e-9 #m (wavelength of light)
P = 13e-6 #W (power of the beam)
waist = 300e-6 #m (beam waist)

#Initial State
F1 = 7/2
mF1 = -7/2
J1 = 1/2
I1 = 4

#Final State
F2 = 9/2
mF2 = -9/2
J2 = 1/2
I2 = 4

#*********************************************************

#Calculate the beam intensity in SI units
I = 2*P/(np.pi*waist**2)

#Calculate the line strength of the electronic transition
S = (3*h*e_0*(2*J2+1)*wavelength**3)/(16*np.pi**3)*gamma
print('Line Strength in a.u.:', S/(e*a_0)**2)

#Use Wigner-Eckart theorem to factor out the mF terms from the dipole matrix elements
WE_term = (-1)**(F2-mF2)*wigner_3j(F2,k,F1,-mF2,q,mF1)

#Use spectator theorem to factor out the nuclear component from the dipole matrix elements
#leaving only the electronic portion (which is given by the line strength of the transition)
if I1 == I2:
    spec_term = (-1)**(J2+I2+F1+k)*np.sqrt((2*F2+1)*(2*F1+1))*wigner_6j(F2,k,F1,J1,I2,J2)
else:
    spec_term = 0
    
#Combine the two terms 
full_term = WE_term*spec_term
print('Coefficient:', full_term)

#Calculate the dipole matrix element
DME = float(full_term*np.sqrt(S))
print('Dipole Matrix Element (C*m):', DME)

#Calculate the expected Rabi frequency
Rabi_freq = DME*np.sqrt(2*I/(c*e_0*hbar**2))
Rabi_MHz = abs(Rabi_freq)/(2*np.pi)/1e6
print('Rabi frequency (MHz):',Rabi_MHz)

