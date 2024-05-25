import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp, arcsin, arccos, arctan2, sinh, cosh
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar, grid
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi




def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#1 Wasserwert


m_kalt_1 = 46.28/1000
m_warm_1 = 57.10/1000

W_wasserwert_1 = matrix("""
0    20.6;
30   20.6;
60   20.6;
90   20.6;
120  20.6;
130  35;
140  34.9;
150  34.7;
160  34.6;
170  34.6;
200  34.4;
230  34.2;
260  34;
290  34
""")# zeit in s, Temp. in gradC

T_h_1 = 48.5


m_kalt_2 = 52.83/1000
m_warm_2 = 54.61/1000

W_wasserwert_2 = matrix("""
0    21.2;
30   21.3;
60   21.4;
90   21.4;
120  21.4;
150  21.4;
180  21.4;
190  33;
200  33;
210  33;
220  32.9;
230  32.8;
240  32.8;
250  32.7;
280  32.6;
210  32.5;
240  32.4;
270  32.3
""")


T_h_2 = 47.7


#2 eis
m_warm_1 = 66.86/1000
m_eis_1 = 17.63/1000
W_eis_1 = matrix("""
0    46.1;
30   45.9;
60   45.5;
90   45.2;
120  44.9;
150  44.7;
180  44.4;
210  44.1;
240  44;
250  27.3;
260  24.3;
270  24.4;
280  24.4;
290  24.4;
300  24.4;
310  24.5;
340  24.6;
370  24.6;
400  24.6;
430  24.7;
460  24.7;
490  24.7
""")



m_warm_2 = 65.36/1000
m_eis_2 = 27.27/1000
W_eis_2 = matrix("""
0    49;
30   48.9;
60   47.4;
90   47;
120  46.5;
150  46.0;
180  45.7;
210  45.4;
240  45;
250  18;
260  16;
270  16;
280  16.2;
290  16.3;
300  16.4;
310  16.5;
340  16.6;
370  16.7;
400  16.9;
430  17;
460  17;
490  17.1
""")
#3.1 kondensationsmethode

m_davor_1 = 103.15/1000
m_danach_1 = 107.61/1000
W_kondens_1 =matrix("""
0   20.6;
30  20.6;
60  22.4;
70  24.0;
80  27.0;
90  29.7;
100 32.7;
110 34.6;
120 37.0;
130 40.0;
140 43.0;
150 45.0;
180 44.0;
210 34.9;
240 34.7;
270 34.5;
""")

m_davor_2 = 94.52/1000
m_danach_2 = 98.61/1000
W_kondens_2 =matrix("""
0   21.8;
30  22;
60  22;
90  25;
100 26.5;
110 28.9;
120 31.5;
130 34.7;
140 37;
150 39.9;
160 42;
170 45;
180 47.7;
190 49.7;
220 48.8;
250 48.6;
280 48.2
""")

#3.2 verdampfungsmethode


#1. durchgang

U1 = 94
I1 = 2
m1 = 18.68/1000

#2. durchgang

U2 = 84
I2 = 1.8
m2 = 15.03/1000

#3. durchgang

U3 = 70
I3 = 1.5
m3 = 10.62/1000

# 4 minuten

#946 hPa druck
#push