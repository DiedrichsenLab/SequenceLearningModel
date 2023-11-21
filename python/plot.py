import matplotlib.pyplot as plt
import numpy as np

def plot_trial(SIM, T, figsize=(10, 15)):
    X=SIM['X']
    B=SIM['B']
    t=SIM['t']
    numOptions, _, numPresses = X.shape

    fg, ax = plt.subplots(nrows=numPresses,ncols=1,figsize=figsize)

    for i in range(numPresses):
        for j in range(numOptions):
            ax[i].plot(t, X[j, :, i], label=f'{j + 1}')
        ax[i].plot(t, B, 'k')
        ax[i].axvline(x=T.stimTime[i], color='r', linestyle=':')
        ax[i].axvline(x=T.decisionTime[i], color='r')
        ax[i].axvline(x=T.pressTime[i], color='k')
        ax[i].legend()

    return fg, ax

def plot_IPI(IPI,label,figsize=(5,3)):

    num = np.arange(1,len(IPI[0])+1)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    colors = ['blue','red','orange','purple','green','yellow','black','pink','brown','gray','cyan','magenta']

    for i in range(len(IPI)):
        ax.plot(num,IPI[i],linewidth=2,label=label[i],color=colors[i])

    ax.set_xticks(num)

    # Remove the box boundary
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set labels and title
    ax.set_ylabel('Inter-Press Interval')
    ax.legend()


    return fig,ax
