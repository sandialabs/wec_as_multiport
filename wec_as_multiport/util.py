# © 2024 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.

from scipy.constants import golden
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from scipy.optimize import fsolve
import numpy as np


def figsize(wf=1, hf=1, columnwidth=250):
    """Parameters:
      - wf [float]:  width fraction in columnwidth units
      - hf [float]:  height fraction in columnwidth units.
                     Set by default to golden ratio.
      - columnwidth [float]: width of the column in latex. Get this from LaTeX 
                             using \showthe\columnwidth
    Returns:  [fig_width,fig_height]: that should be given to matplotlib
    https://stackoverflow.com/questions/29187618/matplotlib-and-latex-beamer-correct-size/30170343
    """

    hf = hf/golden
    fig_width_pt = columnwidth*wf
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*hf      # height in inches
    return [fig_width, fig_height]

def find_maximum_with_interpolation(x_data, y_data, method='cubic', bounds=None):
    """
    Sandia AI
    Finds the maximum value of a function defined by data points using interpolation.

    Parameters:
    - x_data: array-like, the x-coordinates of the data points.
    - y_data: array-like, the y-coordinates of the data points.
    - method: str, the type of interpolation ('linear', 'quadratic', 'cubic', etc.).
    - bounds: tuple, the bounds for the x-values to search for the maximum (min_x, max_x).

    Returns:
    - max_x: the x-coordinate of the maximum point.
    - max_y: the maximum value of the interpolated function.
    """
    
    # Step 1: Interpolate the data
    interpolator = interp1d(x_data, y_data, kind=method, fill_value="extrapolate")

    # Step 2: Define a function for optimization
    def interpolated_function(x):
        return interpolator(x)

    # Step 3: Find the maximum using optimization
    if bounds is None:
        bounds = (min(x_data), max(x_data))
    
    result = minimize_scalar(lambda x: -interpolated_function(x), bounds=bounds, method='bounded')

    # Get the maximum value and corresponding x
    max_x = result.x
    max_y = -result.fun

    return max_x, max_y
        
        
def find_zero_crossings(x, y):
    """Generated from Google AI: Finds zero crossings in a 1D signal.

    Args:
        x (array-like): The x-coordinates of the data points.
        y (array-like): The y-coordinates of the data points.

    Returns:
        array: The x-coordinates of the interpolated zero crossings.
    """

    # Find indices where the sign of y changes
    sign_changes = np.where(np.diff(np.sign(y)))[0]

    # Interpolate the zero crossings
    zero_crossings = []
    for i in sign_changes:
        x1, x2 = x[i], x[i + 1]
        y1, y2 = y[i], y[i + 1]

        # Use linear interpolation to find the zero crossing
        zero_crossing = x1 - y1 * (x2 - x1) / (y2 - y1)
        zero_crossings.append(zero_crossing)

    return np.array(zero_crossings)

def depth_function(k, h=None):
    if h is None:
        h = np.infty
        D = 1
    else:
        D = (1 + 2*k*h/np.sinh(2*k*h))*np.tanh(k*h)
    return D

def w2k(w, h=None, g=9.81):
    """Radial frequency to wave number"""
    if h is None:
        h = np.infty
    func = lambda k: __dispersion__(k=k, w=w, h=h, g=g)
    x0 = w**2/g
    return fsolve(func, x0=x0)[0]

def k2w(k, h=None, g=9.81):
    """Wave number to radial frequency"""
    if h is None:
        h = np.infty
    func = lambda w: __dispersion__(k=k, w=w, h=h, g=g)
    x0 = np.sqrt(g*k)
    return fsolve(func, x0=x0)[0]

def __dispersion__(k, w, h=None, g=9.81):
    """Dispersion relationship"""
    if h is None:
        h = np.infty
    return w**2 - g * k * np.tanh(k * h)
