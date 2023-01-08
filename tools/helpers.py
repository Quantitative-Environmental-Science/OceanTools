import numpy as np

def get_last_values(box):
    """Returns the final values in all time-varying variables in a box.

    Parameters
    ----------
    box : dict
        A dictionary containing the box variables.
        
    Returns
    -------
    dict : containing the last values of all time-varying variables in the box.
    """
    out = {} 
    for k, v in box.items():
        if isinstance(v, np.ndarray):
            out[k] = v[-1]
            
    return out