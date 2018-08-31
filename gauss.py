'''
This script is the main code which imports LGR and AMR data and 
fits the data to a Gaussian plume, then outputs an estimated emission
rate and plots the results.

Author: Colin Arrowsmith
'''

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import utm
from scipy import stats
import gauss_helper as gh
import gauss_plotter as gp


'''  This function takes the input dictionary as an unput and runs the Gaussian
Plume model regression on specific sections of the dataset'''
def main(inp) :
    # Import the data
    data = gh.load_data(inp["file_path"], interval=None)
    background = np.mean(data.ch4.loc[inp["b_interval"][0] : inp["b_interval"][1]])
    data = data.loc[inp["intervals"][0, 0] : inp["intervals"][-1, -1]]

    # Define some variables, calculate and subtract background
    T = np.median(data.T)                           # Average ambient temperature in deg C
    p = np.median(data.p)                           # Average pressure at receptor height in hPa

    data["ch4"] = data["ch4"] - background


    # Convert lats and lons to UTM
    data["utmx"] = data.apply(lambda row: utm.from_latlon(row["lat"], row["lon"])[0], axis=1) 
    data["utmy"] = data.apply(lambda row: utm.from_latlon(row["lat"], row["lon"])[1], axis=1) 
    

    # Loop through transects
    for k in range(len(inp["intervals"])) :
        # Initialize lists
        concentrations = []    # Each source will have it's own array of concentrations
        rates = []             # Each source will have it's own array of betas
        wind_arr = np.zeros([len(inp["sources"]), 3])

        # Spatially average the transect data
        meas_c, x, y = gh.spatial_avg(data.utmx.values, data.utmy.values, data.ch4.values, grid_spacing=inp["grid_size"])
        meas_m = gh.ppm_to_mass(meas_c, T, p)

        # Create array with indices of which spatial interval the receptor falls into
        # Do this only if we wish to split up the data (i.e. spatial_bins != None)
        if  type(inp["spatial_bins"]) != type(None) :
            if inp["axis"] == 0 :
                inds = np.digitize(x, inp["spatial_bins"])

            elif inp["axis"] == 1 :
                inds = np.digitize(y, inp["spatial_bins"])

        elif type(inp["spatial_bins"]) == type(None) :
            inds = np.zeros(np.shape(meas_m))



        # Loop through sources
        for n in range(len(inp["sources"])) :
            # Convert source coordinates to UTM
            source_x, source_y = utm.from_latlon(inp["sources"][n, 0], inp["sources"][n, 1])[0 : 2] 


            # Make an array containing ones where the plume is and NaN elsewhere
            I = np.where(inds == n*2, 1.0, np.nan) * (meas_m/meas_m)

            # Take plume from this source out of meas_c, x, y grids
            m_plume = meas_m * I
            x_plume = x * I
            y_plume = y * I

            # Get ws, wd for transect (comment out if explicitly setting wd, ws)
            # wd, ws = gh.wind_average(data["wd"].loc[inp["w_interval"][n, 0] : inp["w_interval"][n, 1]].values,
            #                          data["ws"].loc[inp["w_interval"][n, 0] : inp["w_interval"][n, 1]].values)
            
            # Explicitly set wind direction and speed for all sources here
            wd, ws = 300.0, 3.0

            # Compute Betas for this plume
            Beta = gh.beta(source_x, source_y, x_plume, y_plume, 
                           inp["z"], inp["H"][n], inp["stability"], ws, wd)

            # Perform linear regression to compute emission rate
            g_per_sec, bgd, r_value, p_value, std_err = stats.linregress(np.nan_to_num(Beta).flatten(), np.nan_to_num(m_plume).flatten())
            kg_per_hr = (g_per_sec/1000)*60*60


            # Calculate downwind concentrations and convert to ppm
            calc_m = g_per_sec * gh.beta(source_x, source_y, x*(meas_m/meas_m), y*(meas_m/meas_m), 
                                         inp["z"], inp["H"][n], inp["stability"], ws, wd)
            calc_c = gh.mass_to_ppm(calc_m, T, p)
            


            # Append to lists
            concentrations.append(calc_c)
            rates.append(kg_per_hr)


            wind_arr[n, 0], wind_arr[n, 1] = source_x, source_y
            wind_arr[n, 2] = np.radians(270.0 - wd)


            print("Transect ", k+1, ", Source", n+1, " : ", kg_per_hr, " kg per h")

        # Plot the results on a map
        gp.gauss_map_plot(meas_c, concentrations, x, y, axis=0, imagery='satellite', transect=k+1, rate_list=rates, wind=wind_arr)
        
        # Plot the results on a regular grid (no map)
        # gp.gauss_plot(meas_c, concentrations, x, y, axis=0, imagery=None, transect=2, rate_list=rates, wind=wind_arr)



    return

if __name__ == "__main__" :

    # Create an input dictionary with all the data needed for the model
    input_dict = {
                  # The path to the csv containing the measured data
                  "file_path": "/home/colin/Documents/GTA-Emissions/Plume-Analysis/Data/sync_data_2018-06-28.csv",

                  # Coordinates of each source
                  "sources": np.array([[43.761588, -79.592489],
                                       [43.762403, -79.583672]]),

                  # The atmospheric stability category assiciated wuth the data
                  "stability": "A",

                  # The date-time intervals for each transect
                  "intervals": np.array([#["2018-06-28 16:01:06", "2018-06-28 16:12:00"],
                                         ["2018-06-28 16:34:00", "2018-06-28 16:45:00"]]),
                                         #["2018-06-28 16:50:00", "2018-06-28 17:11:00"],
                                         #["2018-06-28 17:11:00", "2018-06-28 17:17:00"]]),

                  # Date-time interval to average for wind (one interval for each source).
                  "w_interval": np.array([["2018-06-28 16:52:00", "2018-06-28 16:59:50"],  # Further East
                                          ["2018-06-28 17:02:00", "2018-06-28 17:07:00"]]),# Further West

                  # Time interval in which to calculate the background
                  "b_interval": np.array(["2018-06-28 15:56:08", "2018-06-28 16:01:06"]),

                  # Defines whether measurements primarily transected along Eastings (0) or along Northings (1).
                  "axis": 0,

                  # Spatial subset of data to run seperate regressions on (Whether this is Easting or Northing
                  # is given by "axis")
                  "spatial_bins": np.array([613743.0, 613920.0, 614536.0]),

                  # Effetive plume height of each source above ground (same order as "sources")
                  "H": np.array([3.0, 3.0]),

                  # Receptor height above ground in m (float)
                  "z": 2.0,

                  # Size of averaging grid boxes in m (float)
                  "grid_size": 20.0

                }

    main(input_dict)