# rust-python interface file 
import numpy as np

def date_to_julian_day_num(
    year: int,
    month: int,
    day: int
):
    """
    Gregorian date to julian day

    Inputs
    ------
    year

    month

    day
    """

def calc_earth_day_length(
    lattitude_deg: np.ndarray,
    longitude_deg: np.ndarray,
    julian_day: np.ndarray
) -> np.ndarray :
    """
    Calculate Day length on earth

    Inputs
    ------
    lattitude_deg

    longitude_deg
    
    julian_day
    """
