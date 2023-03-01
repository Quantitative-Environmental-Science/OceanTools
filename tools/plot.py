import matplotlib.pyplot as plt

def boxes(time, vars, *boxes, axs=None, label=None, **kwargs):
    """Plot a set of variables in a set of boxes.

    Parameters
    ----------
    time : array-like
        An array containing the time axis of the model.
    vars : list of str
        A list containing the names of the variables to plot.
    *boxes : dicts
        The boxes that you want to plot. Each box is a dictionary.
    axs : list of Axes, optional
        A list of axes to plot the variables in. If not given, a new figure and axes are created.
    label : str, optional
        Some custom text to add to the legend.
    **kwargs
        Keyword arguments to pass to the pyplot.plot function.
    """
    if isinstance(vars, str):
        vars = [vars]
        
    if axs is None:
        fig, axs = plt.subplots(len(vars), 1, figsize=(8, 1.7 * len(vars)), sharex=True, constrained_layout=True)
    else:
        fig = axs[0].figure

    if hasattr(axs, 'plot'):
        axs = [axs]
        
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
                try:
                    ax.plot(time, box[var], color=cdict[boxname], **kwargs)
                except ValueError:
                    continue
            
        ax.set_ylabel(var)
    
    # box labels, if they're not already there
    current_labels = axs[-1].get_legend_handles_labels()[1]
    plot_orig = False
    for box in boxes:
        if box['name'] not in current_labels:
            axs[-1].plot([], [], color=cdict[box['name']], label=box['name'])
        else:
            plot_orig = True
    axs[-1].legend(fontsize=8)

    if label is not None:
        current_labels = axs[-2].get_legend_handles_labels()[1]
        if plot_orig and 'original' not in current_labels:
            axs[-2].plot([], [], color=(.3,.3,.3), label='original')
        axs[-2].plot([], [], color=(.3,.3,.3), label=label, **kwargs)
        axs[-2].legend(fontsize=8)

    axs[-1].set_xlabel('time')
    axs[-1].set_xlim(0, max(time))
    
    return fig, axs

