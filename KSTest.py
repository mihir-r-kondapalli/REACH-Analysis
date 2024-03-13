# Kolmogorov-Smirnov Test
# Mihir Kondapalli

# Instructions

# ks_test(data1, data2, 0.05, graph = True, print_params = True)
# this is how the test is run

# first dataset, second dataset, significance level

# if graph is true, the test results will be visualized
# default is false

# if print_params is true, the test results will be printed
# default is false

# run matplotlib.pyplot.show() to visualize the results

import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
import scipy.stats as stats
from scipy.stats import ks_2samp as kt

ALPHAS = [0.1, 0.05, 0.025, 0.01, 0.005, 0.001]
COEFFS = [1.22, 1.36, 1.48, 1.63, 1.73, 1.95]

def ks_test(data1, data2, alpha, graph = False, print_params = False, xlabel = 'Data Values', title = 'Kolmorogorov-Smirnov Test'):
    fmt1 = np.sort(np.array(data1))
    fmt2 = np.sort(np.array(data2))
    xf1 = np.unique(fmt1)
    xf2 = np.unique(fmt2)
    cf1 = get_cum_freq(data1)
    cf2 = get_cum_freq(data2)
    
    D, di1, di2 = get_max_dist_approx(xf1, cf1, xf2, cf2)
    
    y2 = 0
    if cf1[di1] > cf2[di2]:
        y2 = cf1[di1]-D
    else:
        y2 = cf1[di1]+D

    x_values = [xf1[di1], xf1[di1]]
    y_values = [cf1[di1], y2]

    if print_params:
        print('Critical KS-Statistic -> ', get_D_crit(len(data1), len(data2), alpha))
        print('Observed KS-Statistic -> ', D)
        print('Significance Level -> ', alpha)
        print('One Sided P-value -> ', get_pval(D, len(data1), len(data2)))
        print()
        #print('One Sided Module P-value -> ', get_pval_module(D, len(data1), len(data2)))
        #print('One Sided KS_2Samp P-value -> ', kt(data1, data2))
    if graph:
        plt.plot(x_values, y_values, 'rx', linestyle="--")
        plt.step(xf1, cf1, where='post')
        plt.step(xf2, cf2, where='post')
        plt.title(title)
        plt.ylabel('Cumulative Frequency')
        plt.xlabel(xlabel)

    return D, get_D_crit(len(data1), len(data2), alpha), get_pval(D, len(data1), len(data2))

def get_cum_freq(data):
    orgd = np.sort(data)
    freqs = []
    length = len(data)
    i = 0
    while i<len(orgd):
        count = (orgd==orgd[i]).sum()
        freqs.append(count/length)
        i+=count-1
        i+=1
    cumfreqs = np.cumsum(np.array(freqs))
    return np.array(cumfreqs)

def get_max_dist_approx(xf1, cf1, xf2, cf2):
    if (xf1[0] < xf2[0] and xf1[len(xf1)-1] < xf2[0]) and (xf1[len(xf1)-1] > xf2[len(xf2)-1] and xf1[0] > xf2[len(xf2)-1]):
        return -1, -1
    D = 0
    D_index1 = 0
    D_index2 = 0
    for i in range(0, len(xf1)):
        if xf1[i] < xf2[0]:
            diff = 0
        elif xf1[i] > xf2[len(xf2)-1]:
            diff = 0
        else:
            index = get_prev_index(xf2, xf1[i])
            diff = abs(cf1[i]-cf2[index])
        if diff>D:
            D = diff
            D_index1 = i
            D_index2 = index
    return D, D_index1, D_index2

def get_prev_index(data, xval):
    
    if xval < data[0]:
        return 0
    elif xval > data[len(data)-1]:
        return len(data)-1
    for i in range(0, len(data)):
        if data[i] > xval:
            return i-1

    return -1

def get_pval_module(D, m, n):
    return stats.kstwo.sf(D, np.round(m*n/(m+n)))

def get_D_crit(n1, n2, alpha):
    try:
        coeff = COEFFS[ALPHAS.index(alpha)]
    except:
        print('Coefficient not found for alpha: ', alpha)
        exit(1)

    return coeff*np.sqrt((n1+n2)/(n1*n2))

def get_pval(D, n1, n2):
    UPPER_LIM = 1000
    z = D*np.sqrt(n1*n2/(n1+n2))
    value = 0   
    for i in range(1, UPPER_LIM):
        value = value + (-1)**(i-1)*np.exp(-2*(i**2)*(z**2))
        #value = value + np.exp(-(2*i-1)**2*(np.pi**2)/(8*(D**2)))
    return 2*value
