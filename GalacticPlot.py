import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

def eq2gal(ra, dec):
    
    '''
    Transforms equatorial coordinates to galactic ones.
    Then prepares them for matplotlib aitoff projection. 
    '''
    
    eq = SkyCoord(ra, dec, unit=u.deg)
    gal = eq.galactic

    # Minus appears because of “mapping from the inside” issue
    l_gal, b_gal = -gal.l.wrap_at('180d').radian, gal.b.radian
    
    return l_gal, b_gal

def ecl2gal(lon_ecl, lat_ecl):
    
    '''
    Transforms ecliptic coordinates to galactic ones.
    Then prepares them for matplotlib aitoff projection.
    '''
    
    ecl = SkyCoord(lon_ecl, lat_ecl, unit=u.deg, frame='barycentricmeanecliptic')
    gal = ecl.transform_to('galactic')

    # Minus appears because of “mapping from the inside” issue
    l_gal, b_gal = -gal.l.wrap_at('180d').radian, gal.b.radian
    
    return l_gal, b_gal

def galactic_plot(ras, decs, sizes, markers, colors, labels, cmaps, ecl = False, title = 'Galactic Plane', scale = 1):
  plt.figure(figsize=(int(scale*14),int(scale*7)))
  plt.subplot(111, projection='aitoff')
  for i in range(0, len(ras)):
    if not(ecl):
        l_eq_gal, b_eq_gal = eq2gal(ras[i], decs[i])
    else:
        l_eq_gal, b_eq_gal = ecl2gal(ras[i], decs[i])
    plt.scatter(l_eq_gal, b_eq_gal, s=sizes[i], marker=markers[i], label=labels[i], c = colors[i], cmap = cmaps[i])
  # Essential thing is to rename RA axis ticks to transform them to conventional format
  plt.xticks(ticks=np.radians([-150, -120, -90, -60, -30, 0, \
                              30, 60, 90, 120, 150]),
            labels=['150°', '120°', '90°', '60°', '30°', '0°', \
                    '330°', '300°', '270°', '240°', '210°'])
  plt.grid(True)
  plt.legend(fontsize=int(scale*16), loc='lower center')
  plt.title(title, fontsize=str(int(scale*16)))
  ax = plt.gca()
  ax.set_facecolor('black')
