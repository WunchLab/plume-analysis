'''
This file contains various helper functions for the Gaussian 
plume analysis python code.

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


'''
This function takes a stability category (one of A, B, C, or D)
as a string and returns the values of the variables a, b, c, d
'''
def stability(category) :
    # Values for x < 1 km
    stability_dict1 = {
        "A": [213.0, 0.894, 440.8, 1.941, 9.270],
        "B": [156.0, 0.894, 106.6, 1.149, 3.300],
        "C": [104.0, 0.894, 61.00, 0.911, 0.000],
        "D": [68.00, 0.894, 33.20, 0.725, -1.70],
        "E": [50.50, 0.894, 22.80, 0.678, -1.30],
        "F": [34.00, 0.894, 14.35, 0.740, -0.35]
    }

    # Values for x > 1 km
    stability_dict2 = {
        "A": [213.0, 0.894, 459.7, 2.094, -9.60],
        "B": [156.0, 0.894, 108.2, 1.098, 2.000],
        "C": [104.0, 0.894, 61.00, 0.911, 0.000],
        "D": [68.00, 0.894, 44.50, 0.516, -13.0],
        "E": [50.50, 0.894, 55.40, 0.305, -34.0],
        "F": [34.00, 0.894, 62.60, 0.180, -48.6]
    }

    for i in stability_dict1 :
        if i == category :
            a1 = stability_dict1[i][0]
            b1 = stability_dict1[i][1]
            c1 = stability_dict1[i][2]
            d1 = stability_dict1[i][3]
            f1 = stability_dict1[i][4]

    for i in stability_dict2 :
        if i == category :
            a2 = stability_dict2[i][0]
            b2 = stability_dict2[i][1]
            c2 = stability_dict2[i][2]
            d2 = stability_dict2[i][3]
            f2 = stability_dict2[i][4]

    return a1, b1, c1, d1, f1, a2, b2, c2, d2, f2


'''
Import the data into a Pandas DataFrame. Takes the data's file path as input.
interval is a tuple, list, or array containing the start and end time of the 
data to use (None by default, meaning all data from that file will be used).  
'''
def load_data(file_path, interval=None) :
    # Define a dataframe from the pertinent data
    df = pd.read_csv(file_path,
                 header=0, names=["gps_time","alt","amb","ch4","ch4d","co",
                                  "co2","co2d","cod","cog","gasp","h2o","hdop",
                                  "heading","lat","lon","pressure","rd1","rd2",
                                  "sog","t","temp","wd_corr","wd_uncorr",
                                  "ws_corr","ws_uncorr"],
                 sep=",",
                 usecols=["gps_time", "lat", "lon", "wd_corr", "ws_corr", "amb", "pressure", 
                          "ch4d", "co2d", "cod", "h2o"],
                 index_col=False, 
                 na_values=("nan", " nan"),
                 na_filter=False)
    

    # Convert date-time string to Pandas DataTimeIndex
    df = df.set_index(pd.DatetimeIndex(df["gps_time"]))

    # Drop NaNs from the dataframe
    df = df.dropna(axis=0)
    df = df.drop("gps_time", axis="columns")
    # Restrict df to the relevant interval
    if interval :
        df = df.loc[interval[0] : interval[1]]

    # Rename some columns
    df = df.rename(columns={"ch4d":"ch4", "co2d":"co2", "cod":"co", 
                            "amb":"T", "pressure":"p", "wd_corr":"wd",
                            "ws_corr":"ws"})

    return df


'''
This function takes an array of wind directions (in degrees clockwise from North)
and an array of wind speeds are outputs an average direction and speed.
'''
def wind_average(wind_d_array, wind_s_array):
    """this function splits an array of measured wind vectors into components
       and outputs their average. The intial wind vectors are relative to
       north as is the functions output
    """
    thetas = np.radians(wind_d_array)
    vs = wind_s_array

    vxs = vs * np.sin(thetas)
    vys = vs * np.cos(thetas)

    vxs_avg = np.average(vxs)
    vys_avg = np.average(vys)

    thetas_avg = np.degrees(np.arctan2(vys_avg, vxs_avg))
    thetas_avg = (450 - thetas_avg) % 360

    vs_avg = np.average(vs)

    return (thetas_avg, vs_avg)




def meter_dist(utm1x, utm1y, utm2x, utm2y, wind_d):
    """accepts two tuples of form (latitude, longitude) and returns 
       distance in x,y coordinates where x direction is downwind of
       source, y direction is perpendicular to the centreline of plume;
       level with the ground."""
    easting_dist = utm2x - utm1x
    northing_dist = utm2y - utm1y
   
    total_dist = np.sqrt((easting_dist ** 2) + (northing_dist ** 2))
    
    theta =  270.0 - wind_d   # Convert angle to an angle counterclockwise from East
    psi = np.degrees(np.arctan2(northing_dist, easting_dist))
    phi = np.radians(theta - psi)
    
    x = total_dist * np.cos(phi)
    y = total_dist * np.sin(phi)

    return (x, y)


'''
This function takes an array of Eastings (utmx) and array of Northings (utmy),
an array of data values to be averaged, then averages the data spatially 
given in grid squares given by the grid_spacing parameter in m. A meshgrid
of spatially averaged values is returned.
'''
def spatial_avg(utmx, utmy, data, grid_spacing=10.0) :
    # Find the bounds of the grid, defined by max and min values of utms
    x_min, x_max = np.nanmin(utmx), np.nanmax(utmx)
    y_min, y_max = np.nanmin(utmy), np.nanmax(utmy)

    # Make regularly spaced x and y utm arrays
    x_arr = np.arange(x_min, x_max+grid_spacing, grid_spacing)
    y_arr = np.arange(y_min, y_max+grid_spacing, grid_spacing)

    # Make x and y meshgrid
    x, y = np.meshgrid(x_arr, y_arr, sparse=False)

    # hist, xedges, yedges = np.histogram2d(utmx, utmy, bins=(x_arr, y_arr))
    H = stats.binned_statistic_2d(utmx, utmy, data, 
                                  statistic="median", 
                                  bins=[x_arr, y_arr])[0]
    return H.T, x[:-1, :-1], y[:-1, :-1]



''' This is the main Gaussian function for the model. Takes position of 
    source and receptor in UTM coordinates, stability, heights, and wind
    and returns the concentration at the receptor location. '''
def beta(sourcex, sourcey, recx, recy, z, H, stability_class, u, wd) :
    sep_vector = meter_dist(sourcex, sourcey, recx, recy, wd)
    x = sep_vector[0]
    y = sep_vector[1]

    # If the source and receptor coordinates are arrays, treat them as such
    if type(x) == float :
        # if receptor is less than 1km downwind, use these values of a,b,c,d,f
        if x < 1000.0 :
            a, b, c, d, f = stability(stability_class)[0 : 5]    
        elif x > 1000.0 :
            a, b, c, d, f = stability(stability_class)[5 : ]

        # Convert x to km and calculate sigmas (these will be in m)
        sigma_y = a * ((x / 1000.0)**b)      # x/1000 has units of km, but sigma has units of m
        sigma_z = c * ((x / 1000.0)**d) + f

    else :
        # Define values for sigma eqns for x > 1km and x < 1km
        a1, b1, c1, d1, f1 = stability(stability_class)[0 : 5]    
        a2, b2, c2, d2, f2 = stability(stability_class)[5 : ]

        # Create ndarray of sigmas so that beta can be calculated without iteration
        sigma_y = np.where(x < 1000.0, a1*((x / 1000.0)**b1), a2*((x / 1000.0)**b2))
        sigma_z = np.where(x < 1000.0, c1*((x / 1000.0)**d1) + f1, c2*((x / 1000.0)**d2) + f2)


    B =  (1.0 / (2.0 * np.pi * sigma_y * sigma_z * u)) * np.exp(-0.5*(y/sigma_y)**2.0) \
          * (np.exp(-0.5*((z - H)/sigma_z)**2.0) + np.exp(-0.5*((z + H)/sigma_z)**2.0))

    return B



def ppm_to_mass(C, temp_deg_C, press_hPa) :
    # function takes ppm CH4 and return g/m^3 CH4
    T = temp_deg_C + 273.15                                     # convert deg C to Kelvin
    N = 6.022e+23                                               # molecues per mol
    m = 16.04                                                   # g/mol CH4
    R = 0.08314                                                 # gas constant (m^3 * hPa * K^-1 * mol^-1)
    mol = 1.66e-18                                              # moles of gas in 10^6 molecules of air
    mass_per_vol = C * ((m / N) * (press_hPa / (mol * R * T)))  # g/m^3
    return mass_per_vol
 
def mass_to_ppm(C, temp_deg_C, press_hPa) :
    # function takes concentration in g/m^3 and returns ppm
    T = temp_deg_C + 273.15                                     # convert deg C to Kelvin
    N = 6.022e+23                                               # molecues per mol
    m = 16.04                                                   # g/mol CH4
    R = 0.08314                                                 # gas constant (m^3 * hPa * K^-1 * mol^-1)
    mol = 1.66e-18 
    ppm = C * ((N/m) * ((mol*R*T)/press_hPa))      
    
    return ppm


''' This function saves the input values of a Gaussian analysis
    as well as the results to one output file. '''
def save_results(input_dict, wd, ws) :
    file = open(file_path, "a+")

    file.write("Data from "+str(input_dict["intervals"][0,0])+
               " to "+str(str(input_dict["intervals"][-1,-1])+
                "\r\n"))
    
    file.write("INPUT PARAMETERS:\r\n")

    file.write("Source Locations:\r\n")
    for i in range(len(input_dict["sources"])) :
        src_num = str(i + 1)
        src_lat, src_lon = str(input_dict["sources"][i, 0]), str(input_dict["sources"][i, 1])
        file.write("Source " + src_num + ": " + src_lat + ", " + src_lon + "\r\n")


    file.write("Stability Class: " + str(input_dict["stability_class"]) + "\r\n")

    file.write("Wind Direction : " + str(wd) + r"$^o$" + "\r\n")
    file.write("Wind Speed: " + str(ws) + " m/s" + "\r\n")

    file.write("Effective Stack Height: " + input_dict["H"] + " m" + "\r\n")

    file.write("Receptor height :" + input_dict["z"])

    


    file.close()
    return

if __name__ == "__main__" :

    # Test spatial averaging function on real data
    path = "/home/colin/Documents/GTA-Emissions/Plume-Analysis/Data/sync_data_2018-06-28.csv"
    data = load_data(path, interval=["2018-06-28 16:01:06", " 2018-06-28 17:20:15"])

    print(data)
    utmx = np.zeros(len(data.lat))
    utmy = np.zeros(len(data.lat))

    for i in range(len(data.lat)) :
        utmx[i], utmy[i] = utm.from_latlon(data.lat[i], data.lon[i])[0 : 2]

    avg, x, y = spatial_avg(utmx, utmy, data.ch4, grid_spacing=50.0)  # 10.0m grid spacing
    plt.figure()
    plt.imshow(avg, extent=[x[0, 0], x[0, -1], y[0, 0], y[-1, 0]], interpolation='None', origin="lower", aspect="equal")
    plt.show()