import numpy as np


class Modifier:
    def __init__(self,
                 boxes,
                 variables,
                 coeff=1,
                 multiply=True
                 ):
        self.variables = variables
        self.boxes = boxes
        self.coeff = coeff
        self.multiply = multiply

    def apply(self, d):
        if d['name'] in self.boxes:
            for var in self.variables:
                try:
                    if self.multiply:
                        d[var] *= self.coeff
                    else:
                        d[var] += self.coeff
                except TypeError:
                    print('Incorrect type modified')
        return d


def modify_single_dict(d, modifiers):
    """ Modifies some variables in a dictionary

    Parameters
    ----------
    d : dict
        Dictionary to be modified
    modifiers : Modifier or array-like
        Modifier or list of modifiers to be applied to the dictionaries

    Returns
    ----------
    d : dict
        return dictionary with modifers applied
    """

    if type(modifiers) == Modifier:
        d = modifiers.apply(d)
    else:
        for modifier in modifiers:
            d = modifier.apply(d)
    return d


def modify_dicts(dicts, modifiers):
    """ Modifies some variables in a list of dicts (makes a copy)

        Parameters
        ----------
        d : dict or array-like
            Dictionary or list of dictionaries to be modified
        modifiers : Modifier or array-like
            Modifier or list of modifiers to be applied to the dictionaries

        Returns
        ----------
        d : dict or array-like
            returns modified dictionaries as individual dict or list of dicts
        """

    dicts = copy_dicts(dicts)

    if type(dicts) == dict:
        new_dicts = modify_single_dict(dicts, modifiers)

    else:
        new_dicts = []
        dicts = copy_dicts(dicts)
        for d in dicts:
            try:
                d = modify_single_dict(d, modifiers)
            except AttributeError:
                print(f'Could Not Modify Dictionary {d["name"]}, due to Attribute Error')
            new_dicts.append(d)

    return new_dicts


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
