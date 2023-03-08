import os
import cmaps
import matplotlib.pylab as pylab
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap
import socket
import cartopy.crs as ccrs
import cartopy
# Size configuration
size = 20
params = {'legend.fontsize': size,
          'figure.figsize': (10, 10),
         'axes.labelsize': size,
         'axes.titlesize': size,
         'xtick.labelsize': size,
         'ytick.labelsize': size,
         'axes.axisbelow' : True}
pylab.rcParams.update(params)

dir_out_variables=os.getenv('dir_out_variables')

instructions = {'W100M': {'contour': 'SLP', 'shading': 'W100M_ano', 'climatology': '30d_rmean_clim'},
                'SOLD': {'contour': 'SLP', 'shading': 'SOLD_ano', 'climatology': '30d_rmean_clim'},
                'T2M': {'contour': 'SLP', 'shading': 'T2M_ano', 'climatology': '30d_rmean_clim'}
}

regime = os.getenv('regime')
variable = os.getenv('variable')
season=os.getenv('season') #relevant for regime
clim_type=os.getenv('clim_type')
lag_inputs = [int(x) for x in os.getenv('lag_inputs').split(';')] #'-4;-2;0;2;4'
interval_input = int(os.getenv('interval_input')) #2

plot_setup = {'T2M': {'unit': "Temperature anomaly [K]", 'domain': [-20, 40, 30, 70], 'boundary': 4.5, 'step': 0.5, 'ticks': [-4, -2, 0, 2, 4], 'bins': 19, 'delta': 15, 'cmap': LinearSegmentedColormap.from_list("name", cmaps.hotcold_18lev(np.hstack([np.linspace(0, 7 / 15), np.linspace(8 / 15, 1)])))},
              'W100M': {'unit': "Wind speed anomaly [m/s]", 'domain': [-20, 40, 30, 70], 'boundary': 5, 'step': 0.1, 'ticks': [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5], 'delta': 15, 'bins': 100, 'cmap': cmaps.BlueWhiteOrangeRed},
              'SOLD': {'unit': "Solar irradiation anomaly [ratio]", 'domain': [-20, 40, 30, 70], 'boundary': 0.4, 'step': 0.005, 'ticks': [-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4], 'delta': 15, 'bins': 160, 'cmap': cmaps.ViBlGrWhYeOrRe},
              'SLP': {'slp_distance': 2, 'range_plot': np.arange(960, 1080, 2), 'linewidth_plot': 0.8}}
for lag_input in lag_inputs:
  dir_composites = "%s/df_composites"%(dir_out_variables)
  shading = xr.open_dataset("%s/%s_%sano_%s_clim_lag%s_%s"%(dir_composites, regime, variable, clim_type, lag_input, lag_input+interval_input))
  contouring = xr.open_dataset("%s/%s_%sabs_lag%s_%s"%(dir_composites, regime, "SLP", lag_input, lag_input+interval_input))/100
  shading = shading.sel(lev=1)
  contouring = contouring.sel(lev=1)
  
  #Plotting via matplotlib functions
  if os.path.isdir("%s/plots/"%(dir_composites))==False: os.system("mkdir %s/plots/"%(dir_composites))
  fig,ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
  lon_min, lon_max, lat_min, lat_max = plot_setup[variable]['domain']
  #Add mapping features
  ax.coastlines(color='darkgray', zorder=3)
  ax.add_feature(cartopy.feature.BORDERS, color='darkgray', zorder=3)
  ax.gridlines(linewidth=0.3, color='lightgray', zorder=2)
  ax.axis('off')
  
  #Shading
  shade = plt.contourf(
      np.arange(lon_min,lon_max+0.5, 0.5), np.arange(lat_min,lat_max+0.5, 0.5),
      shading[variable].sel(lon=slice(lon_min,lon_max), lat=slice(lat_min,lat_max)),
      cmap=plot_setup[variable]['cmap'],
      levels=np.linspace(-plot_setup[variable]['boundary'], plot_setup[variable]['boundary'], num=plot_setup[variable]['bins']),
      vmin=-plot_setup[variable]['boundary'],
      vmax=plot_setup[variable]['boundary'],
      extend="both", zorder=1
      )
  #Contouring = SLP
  cont = plt.contour(
          np.arange(lon_min,lon_max+0.5, 0.5), np.arange(lat_min,lat_max+0.5, 0.5),
      contouring['MSL'].sel(lon=slice(lon_min,lon_max), lat=slice(lat_min,lat_max)),
      levels=plot_setup['SLP']['range_plot'], 
      colors="k", 
      linewidths=plot_setup['SLP']['linewidth_plot'], zorder=4
      )
  plt.clabel(cont, inline=True, fmt="%5.0f", fontsize=14, colors="k")

  plt.savefig("%s/plots/%s_%s_%s_%s_lag%s_%s.pdf"%(dir_composites, variable, regime, 'dunkelflaute', clim_type, lag_input, lag_input+interval_input), bbox_inches='tight', dpi=600)
  plt.clf()
  #plt.show()
  
  #Plot colorbar
  fig2, ax2 = plt.subplots(figsize=(8, 4))
  plt.colorbar(shade, orientation="horizontal", ticks=plot_setup[variable]['ticks'], ax=ax2).ax.set_xlabel(plot_setup[variable]['unit'])
  ax2.remove()
  plt.savefig("%s/plots/%s_colorbar.pdf"%(dir_composites, variable), bbox_inches='tight', dpi=600)
  plt.clf()
  #plt.show()
  