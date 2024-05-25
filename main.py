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


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)


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
sigma_m_kalt_1 = 0.01/1000
m_warm_1 = 57.10/1000
sigma_m_warm_1 = 0.01/1000

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
sigma_T_h_1 = 0.5

T_k_1 = 20.6
sigma_T_k_1 = 0.1

T_m_1 = 35
sigma_T_m_1 = 0.5

cw = 4.192*1000
sigma_cw = 0

K1 = cw*(m_warm_1*(T_h_1 - T_m_1)/(T_m_1 - T_k_1) - m_kalt_1)
sigma_K1 = gauss("cw*(m_warm_1*(T_h_1 - T_m_1)/(T_m_1 - T_k_1) - m_kalt_1)")


t = toArray(W_wasserwert_1[:, 0])
sigma_t = 2*ones(len(t))

T = toArray(W_wasserwert_1[:, 1])
sigma_T = 0.5*ones(len(T))


errorbar(t, T, sigma_T, sigma_t,'x', label='Temperaturmessung')
xlabel('Zeit in s', fontsize=20)
ylabel('Temperatur in °C', fontsize=20)
legend(fontsize=15)
grid()
plt.tight_layout()
savefig('wasserwert1')
show()


m_kalt_2 = 52.83/1000
m_warm_2 = 54.61/1000
sigma_m_kalt_2 = 0.01/1000
sigma_m_warm_2 = 0.01/1000

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
sigma_T_h_2 = 0.5

T_k_2 = 21.4
sigma_T_k_2 = 0.1

T_m_2 = 33
sigma_T_m_2 = 0.5

cw = 4.192*1000
sigma_cw = 0

K2 = cw*(m_warm_2*(T_h_2 - T_m_2)/(T_m_2 - T_k_2) - m_kalt_2)
sigma_K2 = gauss("cw*(m_warm_2*(T_h_2 - T_m_2)/(T_m_2 - T_k_2) - m_kalt_2)")


t = toArray(W_wasserwert_2[:, 0])
sigma_t = 2*ones(len(t))

T = toArray(W_wasserwert_2[:, 1])
sigma_T = 0.5*ones(len(T))


errorbar(t, T, sigma_T, sigma_t,'x', label='Temperaturmessung')
xlabel('Zeit in s', fontsize=20)
ylabel('Temperatur in °C', fontsize=20)
legend(fontsize=15)
grid()
plt.tight_layout()
savefig('wasserwert2')
show()


K = (K1 + K2)/2
sigma_K = gauss("(K1 + K2)/2")


#2 eis
m_warm_1 = 66.86/1000
sigma_m_warm_1 = 0.01/1000
m_eis_1 = 17.63/1000
sigma_m_eis_1 = 0.01/1000

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

def linear(x, a, b):
    return a*x + b


t = toArray(W_eis_1[:, 0])
sigma_t = 2*ones(len(t))

T = toArray(W_eis_1[:, 1])
sigma_T = 0.5*ones(len(T))

errorbar(t, T, sigma_T, sigma_t,'x', label='Temperaturmessung')
plot(t[8:11], T[8:11])
optimizedParameters1, s = opt.curve_fit(linear, t[0:8], T[0:8])
plot(t[0:17], linear(t[0:17], *optimizedParameters1), label="fit1")
optimizedParameters2, s = opt.curve_fit(linear, t[10:], T[10:])
plot(t[4:], linear(t[4:], *optimizedParameters2), label="fit1")
t_zwickel = 248
plt.vlines(t_zwickel, min(T), max(T), label='vline')
xlabel('Zeit in s', fontsize=20)
ylabel('Temperatur in °C', fontsize=20)
legend(fontsize=13, loc='center left')
grid()
plt.tight_layout()
savefig('eis1')
show()

T_h_eis_1 = linear(t_zwickel, *optimizedParameters1)
sigma_T_h_eis_1 = 0.1
sigma_T_k_eis_1 = 0.1
T_k_eis_1 = linear(t_zwickel, *optimizedParameters2)

schmelz_1 = ((T_h_eis_1 - T_k_eis_1)*(cw*m_warm_1 + K) - cw*m_eis_1*T_k_eis_1)/m_eis_1
sigma_schmelz_1 = gauss("((T_h_eis_1 - T_k_eis_1)*(cw*m_warm_1 + K) - cw*m_eis_1*T_k_eis_1)/m_eis_1")



m_warm_2 = 65.36/1000
sigma_m_warm_2 = 0.01/1000
m_eis_2 = 27.27/1000
sigma_m_eis_2 = 0.01/1000
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


t = toArray(W_eis_2[:, 0])
sigma_t = 2*ones(len(t))

T = toArray(W_eis_2[:, 1])
sigma_T = 0.5*ones(len(T))

errorbar(t, T, sigma_T, sigma_t,'x', label='Temperaturmessung')
plot(t[8:11], T[8:11])
optimizedParameters1, s = opt.curve_fit(linear, t[0:8], T[0:8])
plot(t[0:17], linear(t[0:17], *optimizedParameters1), label="fit1")
optimizedParameters2, s = opt.curve_fit(linear, t[10:], T[10:])
plot(t[4:], linear(t[4:], *optimizedParameters2), label="fit1")
t_zwickel = 248
plt.vlines(t_zwickel, min(T), max(T), label='vline')
xlabel('Zeit in s', fontsize=20)
ylabel('Temperatur in °C', fontsize=20)
legend(fontsize=13, loc='center left')
grid()
plt.tight_layout()
savefig('eis2')
show()


T_h_eis_2 = linear(t_zwickel, *optimizedParameters1)
sigma_T_h_eis_2 = 0.1
sigma_T_k_eis_2 = 0.1
T_k_eis_2 = linear(t_zwickel, *optimizedParameters2)

schmelz_2 = ((T_h_eis_2 - T_k_eis_2)*(cw*m_warm_2 + K) - cw*m_eis_2*T_k_eis_2)/m_eis_2
sigma_schmelz_2 = gauss("((T_h_eis_2 - T_k_eis_2)*(cw*m_warm_2 + K) - cw*m_eis_2*T_k_eis_2)/m_eis_2")



#3.1 kondensationsmethode

m_davor_1 = 103.15/1000
sigma_m_davor_1 = 0.01/1000
m_danach_1 = 107.61/1000
sigma_m_danach_1 = 0.01/1000
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
210 43.9;
240 43.7;
270 43.5
""")




t = toArray(W_kondens_1[:, 0])
sigma_t = 2*ones(len(t))

T = toArray(W_kondens_1[:, 1])
sigma_T = 0.5*ones(len(T))

errorbar(t, T, sigma_T, sigma_t,'x', label='Temperaturmessung')
plot(t[1:12], T[1:12])
optimizedParameters1, s = opt.curve_fit(linear, t[0:2], T[0:2])
plot(t[0:17], linear(t[0:17], *optimizedParameters1), label="fit1")
optimizedParameters2, s = opt.curve_fit(linear, t[11:], T[11:])
plot(t[2:], linear(t[2:], *optimizedParameters2), label="fit1")
t_zwickel = 100
plt.vlines(t_zwickel, min(T), max(T), label='vline')
xlabel('Zeit in s', fontsize=20)
ylabel('Temperatur in °C', fontsize=20)
legend(fontsize=15, fancybox=True, framealpha=0.5)
grid()
plt.tight_layout()
savefig('kondens1')
show()


T_h_kondens_1 = linear(t_zwickel, *optimizedParameters2)
sigma_T_h_kondens_1 = 0.1
sigma_T_k_kondens_1 = 0.1
T_k_kondens_1 = linear(t_zwickel, *optimizedParameters1)

T0 = 98
sigma_T0 = 1

kondens1 = ((T_h_kondens_1 - T_k_kondens_1)*(K + cw*m_davor_1) - (T0 - T_h_kondens_1)*(m_danach_1 - m_davor_1))/(m_danach_1 - m_davor_1)
sigma_kondens1 = gauss("((T_h_kondens_1 - T_k_kondens_1)*(K + cw*m_davor_1) - (T0 - T_h_kondens_1)*(m_danach_1 - m_davor_1))/(m_danach_1 - m_davor_1)")


m_davor_2 = 94.52/1000
sigma_m_davor_2 = 0.01/1000
m_danach_2 = 98.61/1000
sigma_m_danach_2 = 0.01/1000
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

t = toArray(W_kondens_2[:, 0])
sigma_t = 2*ones(len(t))

T = toArray(W_kondens_2[:, 1])
sigma_T = 0.5*ones(len(T))

errorbar(t, T, sigma_T, sigma_t,'x', label='Temperaturmessung')
plot(t[2:14], T[2:14])
optimizedParameters1, s = opt.curve_fit(linear, t[0:3], T[0:3])
plot(t[0:17], linear(t[0:17], *optimizedParameters1), label="fit1")
optimizedParameters2, s = opt.curve_fit(linear, t[15:], T[15:])
plot(t[2:], linear(t[2:], *optimizedParameters2), label="fit1")
t_zwickel = 135
plt.vlines(t_zwickel, min(T), max(T), label='vline')
xlabel('Zeit in s', fontsize=20)
ylabel('Temperatur in °C', fontsize=20)
legend(fontsize=13, fancybox=True, framealpha=0.5, loc='upper left')
grid()
plt.tight_layout()
savefig('kondens2')
show()


T_h_kondens_2 = linear(t_zwickel, *optimizedParameters2)
sigma_T_h_kondens_2 = 0.1
sigma_T_k_kondens_2 = 0.1
T_k_kondens_2 = linear(t_zwickel, *optimizedParameters1)

kondens2 = ((T_h_kondens_2 - T_k_kondens_2)*(K + cw*m_davor_2) - (98 - T_h_kondens_2)*(m_danach_2 - m_davor_2))/(m_danach_2 - m_davor_2)
sigma_kondens2 = gauss("((T_h_kondens_2 - T_k_kondens_2)*(K + cw*m_davor_2) - (T0 - T_h_kondens_2)*(m_danach_2 - m_davor_2))/(m_danach_2 - m_davor_2)")




#3.2 verdampfungsmethode


#1. durchgang

U = 94
II = 2
m = 18.68/1000

t = 4*60
sigma_t = 3
sigma_U = 2
sigma_II = 0.02
sigma_m = 0.01

verdampf1 = U*II*t/m
sigma_verdampf1 = gauss("U*II*t/m")

#2. durchgang

U = 84
II = 1.8
m = 15.03/1000

verdampf2 = U*II*t/m
sigma_verdampf2 = gauss("U*II*t/m")

#3. durchgang

U = 70
II = 1.5
m = 10.62/1000


verdampf3 = U*II*t/m
sigma_verdampf3 = gauss("U*II*t/m")
# 4 minuten

#946 hPa druck
#push

#drucken
latexTable(UC(W_wasserwert_1[:, 0], 's'), UC(W_wasserwert_1[:, 1], '^\\circ C'))
latexTable(UC(W_wasserwert_2[:, 0], 's'), UC(W_wasserwert_2[:, 1], '^\\circ C'))

latexTable(UC(W_eis_1[:, 0], 's'), UC(W_eis_1[:, 1], '^\\circ C'))
latexTable(UC(W_eis_2[:, 0], 's'), UC(W_eis_2[:, 1], '^\\circ C'))

latexTable(UC(W_kondens_1[:, 0], 's'), UC(W_kondens_1[:, 1], '^\\circ C'))
latexTable(UC(W_kondens_2[:, 0], 's'), UC(W_kondens_2[:, 1], '^\\circ C'))
