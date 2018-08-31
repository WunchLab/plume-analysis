# README


## What is this code for?
This code performs a _Gaussian Plume Analysis_ on a set of in-situ measurements taken by the 
Wunch Lab LGR-Airmar mobile lab. The code loads the data in Python, averages it into a spatial
grid, performs a linear regression, computes an estimated emission rate for each source, and 
plots the resulting downwind concentrations.


## About the Gaussian model
Please see the following paper for a detailed description of the plume analysis method used in
this code.
_link to final report_


## How to use this code
### Dependencies
These scripts were written and tested in Python 3.5 and have not been tested with other Python
versions. The main computational script requires the Numpy, Matplotlib, Scipy and Pandas packages. 
Another package called `utm` is required and can be accessed [here](https://pypi.org/project/utm/).
One plotting function requires `mpl_toolkits.basemap` in order to make plots on satellite and
street maps. This is not required for the basic computational module.

### Usage
Save `gauss.py`, `gauss_helper.py`, and `gauss_plotter.py` in the same directory. The modeling 
script in `gauss.py` can run without `gauss_plotter.py` if all calls to the plotting functions
are removed. To run it, simply run `gauss.py` from the command line or from a python console.

### The input dictionary
`gauss.py` contains a dictionary object called `input_dict`. This contains all the input data
that the model needs. Each key in the dictionary is described below:

* `file_path`: (str). This is the relative or absolute path (including file name) to the file containing the raw data. The file must be a csv with headers that have the same naming structure as the data files [here](https://dataverse.scholarsportal.info/dataset.xhtml?persistentId=doi:10.5683/SP/DEQJGQ).

* `sources`: (np.array). This is a 2-D numpy array containing the coordinates of each source. The array must have the form:
        
```python
np.array([[source1_lat, source1_lon],
          [source2_lat, source2_lon]])
```

* `stability`: (str) This is the atmospheric stability category as defined by Turner, 1970. Can be one of A, B, C, D, E, or F.

* `intervals`: (np.array of datatime strings). The datetime interval for each transect. The model can process several transects seperately, so this is where we define the data that make up each transect. The following is an example of this array:

```python
np.array([["2018-06-28 16:01:06", "2018-06-28 16:12:00"],
          ["2018-06-28 16:34:00", "2018-06-28 16:45:00"],
          ["2018-06-28 16:50:00", "2018-06-28 17:11:00"],
          ["2018-06-28 17:11:00", "2018-06-28 17:17:00"]])
```

* `w_interval`: (np.array of datetime strings). The datetime interval in which to calculate the average wind speed and direction for each source. Must have the same length as `sources`.

* `b_interval`: (np.array of datetime strings). The datetime interval in which to calculate the background pollutant concentration. Must have shape (2, ).

* `axis`: (int). Defines axis along which to average for plots. If 0, plots will be averaged along Northings and concentration plotted as a function of Easting. If 1, plots will be averaged along Eastings and concentration plotted as a function of Norting.

* `spatial_bins`: (np.array of floats or `None`). Eastings or Northings along which to split up data for separate regressions. Eastings must be used if axis=0 and Northings must be used if axis=1. This divides the data into spatial bins so that each source can be treated independently. It is assumed that the plume from source 1 is located in values less than the first entry, the plume from source 2 is located between the 2nd and 3rd entries, and so on. For example, if `"spatial_bins": np.array([613743.0, 613920.0, 614536.0])` and `"axis": 0`, then data located to the West of 613743.0 will be treated as the plume from source1, data between 613743.0 and 613920.0 will be assumed to be between the plumes, and data between 613920.0 and 614536.0 will be assumed to be the plume from source2. If there is only one source and one plume, this should be `None`.

* `H`: (np.array of floats). The effective plume height of each source in meters. Must have the form `np.array([source1_H, source2_H, ...])`.

* `grid_size`: (float). The size of spaially averaged grid boxes in meters.
