'''
This script contains several functions which aid in plotting
the results from the gauss.py code. These are useful for plotting
spatial concentration data on a map or grid.

Author: Colin Arrowsmith
'''

# Import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
matplotlib.style.use('ggplot')
import utm

def make_map(xmin, xmax, ymin, ymax, axes, imagery) :

    # Convert from UTM to latlon
    LLcoord = utm.to_latlon(xmin, ymin, 17, 'T')
    URcoord = utm.to_latlon(xmax, ymax, 17, 'T')

    m = map_prep(LLcoord, URcoord, imagery, axes)


    return m




'''
make_map() returns a matplotlib map object.
LLcoord: Tuple, array, or list of coordinates of lower left corner of map. Has the form (lat, lon).
URcoord: Tuple, array, or list of coordinates of upper right corner of map. Has the form (lat, lon).
layer: string specifying the map tile to be used. Either "satellite" or "streets".
'''
def map_prep(LLcoord, URcoord, layer, axes) :
    resolution = 700
    EPSG = 3978
    map1 = Basemap(llcrnrlon=LLcoord[1], llcrnrlat=LLcoord[0], 
                   urcrnrlon=URcoord[1], urcrnrlat=URcoord[0], 
                   epsg=EPSG, ax=axes, anchor="NW"
                   )

    if layer == "satellite" :
        s = 'ESRI_Imagery_World_2D'
    elif layer == "streets" :
        s = 'ESRI_StreetMap_World_2D'
    else :
        return "Not a valid imagery type.\n Imagery type must be 'satellite' or 'streets'."


    map1.arcgisimage(service=s, verbose= True, xpixels=1000, dpi=resolution, ax=axes)
    
    return map1



'''  This function takes the resuts of a Gaussian plume analysis and plots:
1.  The measured concentrations, and modeled concentrations due to each source
    averaged along one axis (easting or northing).
2.  The measured and modelled concentrations on a grid, represented by colored 
    grid boxes.     
- "axis" determines whether to average along Northings (0) or Eastings (1)

- "wind" is an array of the form [[source1x, source1y, source1_wd], [source2x, source2y, source2_wd]...]
'''
def gauss_plot(meas_c, model_c_list, x_grid, y_grid, axis=0, transect=None, rate_list=None, wind=None) :

    # Set up axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 2), (1, 0))
    ax3 = plt.subplot2grid((2, 2), (1, 1))



    # Make plots along axis
    for i in range(len(model_c_list)) :
        if axis == 0:
            dist = x_grid[0, :]
            xlabel = "Easting (m)"
        else :
            dist = y_grid[:, 0]
            xlabel = "Northing (m)"

        # Add emission rate text
        if rate_list :
            fig.text(0, 0-i*0.05+len(rate_list)*0.03, "Source "+str(i+1)+" : "+str(rate_list[i])+" kg/h",
                     fontsize=12)

        # Add arrows to source
        try :
            ax2.quiver(wind[i, 0], wind[i, 1], np.cos(wind[i, 2]), np.sin(wind[i, 2]), scale=10.0, angles="xy", zorder=100)
            ax3.quiver(wind[i, 0], wind[i, 1], np.cos(wind[i, 2]), np.sin(wind[i, 2]), scale=10.0, angles="xy", zorder=100)
        except :
            TypeError

        conc = np.nanmean(model_c_list[i], axis=axis)
        ax1.plot(dist, conc, label="Source "+str(i+1))
    
    ax1.plot(dist, np.nanmean(meas_c, axis=axis), label="Measured")
    ax1.plot(dist, np.nanmean(sum(model_c_list), axis=axis), label="Sum")
    ax1.legend()
    ax1.set_title("Transect 2")#+str(transect))
    ax1.set_ylabel(r"$CH_4$ (ppmv)")
    ax1.set_xlabel(xlabel)

    # Make xy plots of measured data
    ax2.imshow(meas_c, 
               extent=[x_grid[0, 0], x_grid[0, -1], y_grid[0, 0], y_grid[-1, 0]],
               interpolation='None', origin="lower", aspect="equal", zorder=100)
    ax2.set_title("Measured")
    ax2.set_xlabel("Easting (m)")
    ax2.set_ylabel("Northing (m)")

    # Make xy plots of modeled data
    ax3.imshow(sum(model_c_list), 
               extent=[x_grid[0, 0], x_grid[0, -1], y_grid[0, 0], y_grid[-1, 0]],
               interpolation='None', origin="lower", aspect="equal", zorder=100)
    ax3.set_title("Modeled")
    ax3.set_xlabel("Easting (m)")
    ax3.set_ylabel("Northing (m)")


    plt.show()
    return




'''This function plots the same data as gauss_plot, but plots the bottom two figures on 
   a map. imagery can be either 'satellite' or 'streets'. '''
def gauss_map_plot(meas_c, model_c_list, x_grid, y_grid, axis=0, imagery='satellite', transect=None, rate_list=None, wind=None) :

    

    # Convert x, y UTM grid to latlon
    lats = np.zeros(np.shape(x_grid))
    lons = np.zeros(np.shape(y_grid))

    for k in range(len(y_grid)) :
        for j in range(len(x_grid[0, :])) :
            lats[k, j], lons[k, j] = utm.to_latlon(x_grid[k, j], y_grid[k, j], 17, 'T')




    # Determine spatial bounds of the map area
    xarr = x_grid[0, :]
    yarr = y_grid[:, 0]
    try :
        xarr = np.concatenate(xarr, wind[:, 0])
        yarr = np.concatenate(yarr, wind[:, 1])
    except: 
        ValueError
    ymin, ymax = min(yarr), max(yarr)
    xmin, xmax = min(xarr), max(xarr)


    # Set up axes
    fig = plt.figure()
    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 2), (1, 0))
    ax3 = plt.subplot2grid((2, 2), (1, 1))

    # Create map objects
    m1 = make_map(xmin, xmax, ymin, ymax, ax2, imagery)
    m2 = make_map(xmin, xmax, ymin, ymax, ax3, imagery)


    # Make plots along axis
    for i in range(len(model_c_list)) :
        if axis == 0:
            dist = x_grid[0, :]
            xlabel = "Easting (m)"
        else :
            dist = y_grid[:, 0]
            xlabel = "Northing (m)"

        # Add emission rate text
        if rate_list :
            fig.text(0, 0-i*0.05+len(rate_list)*0.03, "Source "+str(i+1)+" : "+str(rate_list[i])+" kg/h",
                     fontsize=12)

        # Add arrows to source
        try :
            windlat, windlon = utm.to_latlon(wind[i, 0], wind[i, 1], 17, 'T')
            windx, windy = m1(windlon, windlat)
            m1.quiver(windx, windy, np.cos(wind[i, 2]), np.sin(wind[i, 2]), scale=10.0, angles="xy", zorder=1)
            m2.quiver(windx, windy, np.cos(wind[i, 2]), np.sin(wind[i, 2]), scale=10.0, angles="xy", zorder=5)
        except :
            TypeError

        conc = np.nanmean(model_c_list[i], axis=axis)
        ax1.plot(dist, conc, label="Source "+str(i+1))


    # Plot the rest of the data on the first figure
    ax1.plot(dist, np.nanmean(meas_c, axis=axis), label="Measured")
    ax1.plot(dist, np.nanmean(sum(model_c_list), axis=axis), label="Sum")
    ax1.legend()
    ax1.set_title("Transect 2")#+str(transect))
    ax1.set_ylabel(r"$CH_4$ (ppmv)")
    ax1.set_xlabel(xlabel)

    # Convert lats and lons to ax coordinates
    x, y = m1(lons, lats)


    # Make xy plots of measured data
    meas_c = np.ma.masked_where(np.isnan(meas_c), meas_c)
    meas = m1.pcolormesh(x, y, meas_c, antialiased=True, ax=ax2)
    ax2.set_title("Measured")
    ax2.set_xlabel("Easting (m)")
    ax2.set_ylabel("Northing (m)")
    m1.colorbar(meas, location='right', label=r"$CH_4$ (ppmv)")

    # Make xy plots of modeled data
    modl_c = np.ma.masked_where(np.isnan(sum(model_c_list)), sum(model_c_list))
    modl = m2.pcolormesh(x, y, modl_c, antialiased=True, ax=ax3)
    m2.colorbar(modl, location="right", label=r"$CH_4$ (ppmv)")
    ax3.set_title("Modeled")
    ax3.set_xlabel("Easting (m)")
    ax3.set_ylabel("Northing (m)")

    # Add arrow pointing North
    ax2.quiver(100, 1010, np.cos(np.radians(102)), np.sin(np.radians(102)), scale=15.0, color="r")
    ax2.annotate('N', xy=(80, 950), xycoords='data', color='r')
    ax3.quiver(100, 1010, np.cos(np.radians(102)), np.sin(np.radians(102)), scale=15.0, color="r")
    ax3.annotate('N', xy=(80, 950), xycoords='data', color='r')
    plt.show()

    return