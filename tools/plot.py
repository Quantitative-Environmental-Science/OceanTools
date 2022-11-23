import matplotlib.pyplot as plt

def boxes(time, vars, *boxes, axs=None, label=None, **kwargs):
    """Plot a set of variables in a set of boxes.

    Parameters
    ----------
    time : array-like
        An array containing the time axis of the model.
    vars : list of str
        A list containing the names of the variables to plot.
    axs : list of Axes, optional
        A list of axes to plot the variables in. If not given, a new figure and axes are created.
    *boxes : dicts
        The boxes that you want to plot. Each box is a dictionary.
    **kwargs
        Keyword arguments to pass to the plot function.
    """
    if axs is None:
        fig, axs = plt.subplots(len(vars), 1, figsize=(8, 1.7 * len(vars)), sharex=True, constrained_layout=True)
    else:
        fig = axs[0].figure

    cdict = {
        'deep': 'tab:grey',
        'hilat': 'tab:blue',
        'lolat': 'tab:red',
        'atmos': 'tab:orange',
    }
    
    for var, ax in zip(vars, axs):
        for box in boxes:            
            if var in box:
                boxname = box['name'] 
                ax.plot(time, box[var], color=cdict[boxname], **kwargs)
            
        ax.set_ylabel(var)

    for box in boxes:
        if label is None:
            line_label = box['name']
        else:
            line_label = box['name'] + ': ' + label
        axs[-1].plot([], [], color=cdict[box['name']], label=line_label, **kwargs)

    axs[-1].legend()
    axs[-1].set_xlabel('time')
    axs[-1].set_xlim(0, max(time))
    
    return fig, axs
