import numpy as np 
from mpl_toolkits.basemap import Basemap 
import matplotlib.pyplot as plt 
import sys 

# from http://www.xsi-blog.com/archives/115

def fibonacci(N): 
    inc = np.pi * (3 - np.sqrt(5)) 
    off = 2. / N 
    r2d = 180./np.pi 
    k = np.arange(0,N) 
    y = k*off - 1. + 0.5*off 
    r = np.sqrt(1 - y*y) 
    phi = k * inc 
    x = np.cos(phi)*r 
    z = np.sin(phi)*r 
    theta = np.arctan2(np.sqrt(x**2+y**2),z) 
    phi = np.arctan2(y,x) 
    lats = 90.-r2d*theta 
    lons = r2d*phi 
    return lats, lons 

npts = int(sys.argv[1]) 

lats, lons = fibonacci(npts) 


print lats,lons


map = Basemap(projection ='ortho',lat_0=0,lon_0=-90) 
map.drawcoastlines() 
map.fillcontinents(color='coral') 
map.drawmapboundary(fill_color='aqua') 
x,y = map(lons, lats) 
map.scatter(x,y,10,marker='o',zorder=10) 
plt.show() 
