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
p     = 4
r_ref = 4.2
D_e   = 6599.300
T_e   = 8832.97
r_e   = 4.179865
a = [
 0.56545546,
 0.11075633,
 0.07116331,
 0.01384038,
-0.13620815,
-0.09931222,
 1.32797985,
 0.86468958,
-3.57895961,
-2.60387730,
 6.53483765,
 2.71676634,
-5.19805450,
-0.06610630,
 0.07355150
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



