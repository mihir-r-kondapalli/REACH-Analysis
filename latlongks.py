from KSTest2D import ks2d2s
from KSTest import ks_test
from GalacticPlot import galactic_plot
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def column(matrix, i):
    return [row[i] for row in matrix]

def run_2samp_test(data1, data2):
    return ks2d2s(column(data1, 1), column(data1, 2), column(data2, 1),  column(data2, 2), extra = True)

def run_2samp_test(data1x, data1y, data2x, data2y):
    return ks2d2s(data1x, data1y, data2x,  data2y, extra = True)

def get_lin_dist(start, stop, length):
    lindist = []
    m = (stop-start)/(length-1)
    for i in range(0, length):
        lindist.append(m*i+start)
    return lindist

f_coords = np.loadtxt('Data Files/FRB_coords.csv', delimiter = ',', skiprows = 1)
m_coords = np.loadtxt('Data Files/Mag_coords.csv', delimiter = ',', skiprows = 1)
s_coords = np.loadtxt('Data Files/Short_coords.csv', delimiter = ',', skiprows = 1)
sv_coords = np.loadtxt('Data Files/ValidShort_coords.csv', delimiter = ',', skiprows = 1)
fv_coords = np.loadtxt('Data Files/ValidFRB_coords.csv', delimiter = ',', skiprows = 1)

sv_lat = column(sv_coords, 1)
sv_long = column(sv_coords, 2)
sv_redshifts = column(sv_coords, 3)
sv_ra = column(sv_coords, 4)
sv_dec = column(sv_coords, 5)

fv_lat = column(fv_coords, 1)
fv_long = column(fv_coords, 2)
fv_redshifts = column(fv_coords, 3)
fv_ra = column(fv_coords, 4)
fv_dec = column(fv_coords, 5)

#D, DC, P = ks_test(column(f_coords, 1), column(m_coords, 1), 0.05, graph = True, print_params = True, xlabel = 'Latitude', title = 'FRBs vs Magnetars')
#plt.show()
#D, DC, P = ks_test(column(f_coords, 2), column(m_coords, 2), 0.05, graph = True, print_params = True, xlabel = 'Longitude', title = 'FRBs vs Magnetars')
#plt.show()
#D, DC, P = ks_test(column(f_coords, 1), column(s_coords, 1), 0.05, graph = True, print_params = True, xlabel = 'Latitude', title = 'FRBs vs Short GRBs')
#plt.show()
#D, DC, P = ks_test(column(f_coords, 2), column(s_coords, 2), 0.05, graph = True, print_params = True, xlabel = 'Longitude', title = 'FRBs vs Short GRBs')
#plt.show()


#p, D = run_2samp_test(f_coords, m_coords)
#print('P-value --> ', p)
#print('D-value --> ', D)
#print()
#p, D = run_2samp_test(f_coords, s_coords)
#print('P-value --> ', p)
#print('D-value --> ', D)
#print()


def get_adj_latlong(compare, lats, longs, ras, decs, limlow, limhigh):
    new_lats = []
    new_longs = []
    new_ras = []
    new_decs = []

    for i in range(0, len(compare)):
        if compare[i]>=limlow and compare[i]<=limhigh:
            new_lats.append(lats[i])
            new_longs.append(longs[i])
            new_ras.append(ras[i])
            new_decs.append(decs[i])

    return new_lats, new_longs, new_ras, new_decs

def print_regions():
    MIN_REDSHIFT = 0
    MAX_REDSHIFT = np.max(sv_redshifts)
    MIN_NUM = 20
    step = 0.05
    bounds_decimals = 2
    p_decimals = 4

    start = MIN_REDSHIFT
    stop = MIN_REDSHIFT+step
    
    while start < MAX_REDSHIFT:
        print()
        print('---'+str(round(start, bounds_decimals))+'---')
        print()
        while stop <= MAX_REDSHIFT:
            new_sv_lats, new_sv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(sv_redshifts, sv_lat, sv_long, sv_ra, sv_dec, start, stop)
            new_fv_lats, new_fv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(fv_redshifts, fv_lat, fv_long, fv_ra, fv_dec, start, stop)
            length_s = len(new_sv_lats)
            length_f = len(new_fv_lats)
            if length_s >= MIN_NUM and length_f >= MIN_NUM:
                D1, DC1, P1 = ks_test(new_fv_lats, new_sv_lats, 0.05, graph = False, print_params = False)
                D2, DC2, P2 = ks_test(new_fv_longs, new_sv_longs, 0.05, graph = False, print_params = False)
                print('['+str(round(start, bounds_decimals))+','+str(round(stop, bounds_decimals))+'] -> ('
                      + str(round(P1, p_decimals))+' , '+str(round(P2, p_decimals))+') : <'+str(length_f)+','+str(length_s)+'>')
            stop+=step
        start+=step
        stop = start + step

#print_regions()

# Prints the region with the most uniform distribution of FRBs

def print_max_uniformity():
    MIN_REDSHIFT = 0
    MAX_REDSHIFT = np.max(sv_redshifts)
    MIN_NUM = 20
    step = 0.02
    bounds_decimals = 2
    p_decimals = 4

    start = MIN_REDSHIFT
    stop = MIN_REDSHIFT+step

    P_lat_max = 0
    P_long_max = 0
    P_max_start = -1
    P_max_stop = -1
  
    while start < MAX_REDSHIFT:
        while stop <= MAX_REDSHIFT:
            new_sv_lats, new_sv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(sv_redshifts, sv_lat, sv_long, sv_ra, sv_dec, start, stop)
            new_fv_lats, new_fv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(fv_redshifts, fv_lat, fv_long, fv_ra, fv_dec, start, stop)
            length_s = len(new_sv_lats)
            length_f = len(new_fv_lats)
            if length_s >= MIN_NUM and length_f >= MIN_NUM:
                D1, DC1, P1 = ks_test(new_fv_lats, get_lin_dist(np.min(new_fv_lats), np.max(new_fv_lats), 100), 0.05, graph = False, print_params = False)
                D2, DC2, P2 = ks_test(new_fv_longs, get_lin_dist(np.min(new_fv_longs), np.max(new_fv_longs), 100), 0.05, graph = False, print_params = False)
                Psum = P1 + P2
                if Psum > (P_lat_max+P_long_max):
                    P_lat_max = P1
                    P_long_max = P2
                    P_max_start = start
                    P_max_stop = stop
            stop+=step
        start+=step
        stop = start + step

    print('Maximum Uniformity --> ['+str(round(P_max_start, bounds_decimals))+','+str(round(P_max_stop, bounds_decimals))+'] : <'
          +str(round(P_lat_max, p_decimals))+','+str(round(P_long_max, p_decimals))+'>')

    print('\n\n')
    print('------------------')
    print('\n\n')

    return P_max_start, P_max_stop

def get_bias(start):
    MIN_REDSHIFT = start
    MAX_REDSHIFT = np.max(sv_redshifts)
    MIN_NUM = 20
    step = 0.02
    bounds_decimals = 2
    p_decimals = 4

    start = MIN_REDSHIFT
    stop = MIN_REDSHIFT+step

    P_lat_max = 0
    P_long_max = 0
    P_max_start = -1
    P_max_stop = -1

    stops = []
    ps = []
    nums = []
    
    while stop <= MAX_REDSHIFT:
        new_sv_lats, new_sv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(sv_redshifts, sv_lat, sv_long, sv_ra, sv_dec, start, stop)
        new_fv_lats, new_fv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(fv_redshifts, fv_lat, fv_long, fv_ra, fv_dec, start, stop)
        length_s = len(new_sv_lats)
        length_f = len(new_fv_lats)
        if length_s >= MIN_NUM and length_f >= MIN_NUM:
            D1, DC1, P1 = ks_test(new_fv_lats, get_lin_dist(np.min(new_fv_lats), np.max(new_fv_lats), 100), 0.05, graph = False, print_params = False)
            D2, DC2, P2 = ks_test(new_fv_longs, get_lin_dist(np.min(new_fv_longs), np.max(new_fv_longs), 100), 0.05, graph = False, print_params = False)
            Ptot = P1 * P2
            stops.append(stop)
            ps.append(1-Ptot)
            nums.append(length_f)

            print('[0, '+str(stop)+'] --> '+str(Ptot))
            
        stop+=step

    return stops, ps, nums

#pstart, pstop = print_max_uniformity()


stops, ps, nums = get_bias(0)
plt.step(stops, ps, where='post')
plt.title('Observational Bias vs Redshift of FRBs')
plt.xlabel('Maximum Redshift')
plt.ylabel('Observational Bias')
plt.show()

'''plt.step(stops, nums, where='post')
plt.title('Number of FRBs vs Redshift of FRBs')
plt.xlabel('Maximum Redshift')
plt.ylabel('Number of FRBs')

plt.show()'''

'''pstart, pstop = 0.55, 2.55

new_fv_lats, new_fv_longs, new_fv_ras, new_fv_decs = get_adj_latlong(fv_redshifts, fv_lat, fv_long, fv_ra, fv_dec, pstart, pstop)
new_sv_lats, new_sv_longs, new_sv_ras, new_sv_decs = get_adj_latlong(sv_redshifts, sv_lat, sv_long, sv_ra, sv_dec, pstart, pstop)

print('2D KS Test:')
p, D = run_2samp_test(new_fv_lats, new_fv_longs, new_sv_lats, new_sv_longs)
print('P-value --> ', p)
print('D-value --> ', D)
print()

galactic_plot([new_fv_ras, new_sv_ras], [new_fv_decs, new_sv_decs], [8, 20], ['v', 'o'],
              ['yellow', 'green'], ['FRBs', 'Short GRBs'], [None, None], ecl = False,
              title = 'Galactic Plane: ['+str(round(pstart, 2))+','+str(round(pstop, 2))+']',
              scale = 0.75)
plt.show()

D, DC, P = ks_test(new_fv_lats, new_sv_lats, 0.05, graph = True, print_params = True, xlabel = 'Latitude', title = 'FRBs vs Short GRBs')
plt.show()
D, DC, P = ks_test(new_fv_longs, new_sv_longs, 0.05, graph = True, print_params = True, xlabel = 'Longitude', title = 'FRBs vs Short GRBs')
plt.show()'''
