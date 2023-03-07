import numpy as np

def get_last_values(*boxes):
    """Returns the final values in all time-varying variables in a box.

    Parameters
    ----------
    box : dict(s)
        One or more model boxes.
        
    Returns
    -------
    dict : If a single box is given, returns a dictionary containing the final values of all time-varying attributes. If multiple boxes are given, returns a dictionary containing the final values of all time-varying attributes for each box.
    """
    out = {} 
    
    if len(boxes) == 1:
        for k, v in boxes[0].items():
            if isinstance(v, np.ndarray):
                out[k] = v[-1]        
        return out
    else:
        for box in boxes:
            out[box['name']] = get_last_values(box)
    
    return out


def copy_dicts(dicts):
    """Returns a copy of each dictionary in the list (or of an individual dictionary)

    Parameters
    ----------
    dicts : array-like or dict
        array of dictionaries to be copied

    Returns
    -------
    new_dicts : list
        list of copied dictionaries
    """

    if type(dicts) == dict:
        return dicts.copy()

    else:
        try:
            new_dicts = []
            for dictionary in dicts:
                new_dicts.append(dictionary.copy())

            return new_dicts

        except AttributeError:
            print('Dictionaries could not be copied, output was the original dicts')
            return dicts
