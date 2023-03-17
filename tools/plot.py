import matplotlib.pyplot as plt


def boxes(time, vars, *boxes, axs=None, label=None, height=1.7, width=8, **kwargs):
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
    height : float, optional
        Height modifier to edit the height of the figure
    width : float, optional
        Width modifier to edit the width of the figure
    **kwargs
        Keyword arguments to pass to the pyplot.plot function.
    """
    if isinstance(vars, str):
        vars = [vars]

    if axs is None:
        fig, axs = plt.subplots(len(vars), 1, figsize=(width, (height * len(vars))), sharex=True,
                                constrained_layout=True)
    else:
        fig = axs[0].figure
        
    if hasattr(axs, 'plot'):
        axs = [axs]

    if not hasattr(fig, 'n_layers'):
        fig.n_layers = 0
        fig.box_labels = []
        fig.line_labels = []
    
    # number of panels
    n_plots = len(axs)

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
    if fig.n_layers == 0:
        for box in boxes:
            line = axs[-1].plot([], [], color=cdict[box['name']], label=box['name'])
            fig.box_labels += line
    else:
        current_labels = axs[-1].get_legend_handles_labels()[1]
        if box['name'] not in current_labels:
            line = axs[-1].plot([], [], color=cdict[box['name']], label=box['name'])
            fig.box_labels += line
    axs[-1].legend(fontsize=8)

    # line labels
    if n_plots == 1:
        label_ax = axs[-1]
    else:
        label_ax = axs[-2]

    if (fig.n_layers > 0) and (len(fig.line_labels) == 0):
        line = label_ax.plot([], [], color=(.3, .3, .3), label='original')
        fig.line_labels += line
        
    if label is not None:
        line = label_ax.plot([], [], color=(.3, .3, .3), label=label, **kwargs)
        fig.line_labels += line
        
        label_ax.legend(fontsize=8)
    
    axs[-1].set_xlabel('time')
    axs[-1].set_xlim(0, max(time))

    fig.n_layers += 1
    
    return fig, axs
