########################################################################
#
# Implementation of the Expanded Morse Oscillator (EMO) potential
#
# E. G. Lee, J. Y. Seto, T. Hirao, P. F. Bernath, R. J. LeRoy,
# J. Mol. Spectrosc. 194, 197 (1999).
#
########################################################################
import sys
import math

########################################################################
#
#                   PARAMETERS OF THE POTENTIAL
#
########################################################################
# M. Tamanis, I. Klincare, A. Kruzins, O. Nikolayeva, R. Ferber,
# E. A. Pazyuk and A. V. Stolyarov.
# Direct excitation of the "dark" b3Pi state predicted by deperturbation
# analysis of the A1Sigma+ -- b3Pi complex in KCs.
# Phys. Rev. A 82, 032506 (2010)
########################################################################
r_min = 3.5
r_max = 7.0
step  = 0.1
p     = 3
r_ref = 5.0
D_e   = 5567.546
T_e   = 10049.40
r_e   = 4.981380
a = [
 0.44694080,
 0.01225783,
 0.01194530,
 0.11642869,
 0.16605806,
 0.40208001,
-0.27784050,
-1.11793578,
 1.10616209,
 1.82170820,
-1.21749951,
-1.14724244,
 0.00574038
]
########################################################################

def pot_emo(r):
	y = (r**p - r_ref**p) / (r**p + r_ref**p)
	alpha = 0.0
	for i,a_i in enumerate(a):
		alpha += a_i * (y**i)
	return T_e + D_e*(1-math.exp(-alpha*(r-r_e)))**2

r = r_min
while r <= r_max:
	print('%16.8f%30.16f' % (r, pot_emo(r)))
	r += step



