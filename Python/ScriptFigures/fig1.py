import numpy as np
import fmatoolbox as fma
import isautilities as isau
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig1'
LEFT = 0.7; TOP = 0.3
regs = ('nr','hpc','pfc')
states = ('sws','rem','other')
state_names = ('nREM','REM','wake')
fig = fma.plotting.figure(n=0, size=(17.5,15))
results = isau.load_results(res_root / 'exampleISA.pkl')


# example raster
ax = isau.axes(fig, LEFT, TOP+0.5,7.3,5)
start = results['start']
stop = results['stop']
label = ('firing rate (Hz)', None, None)
ticks = ([0,1,2], [], [])
n = [[0,0]]
cbs = []
for i, r in enumerate(regs):
    fr = results[f'fr_{r}']
    n_units = results[f'units_{r}']
    t = fr[:,0]
    fr = fr[:,1:]
    n.append([n[-1][1], n[-1][1] + n_units])
    im, cb = fma.plotting.plotColorMap(fr.T, x=t, y=n[-1], vmin=0, vmax=2, cmap=LinearSegmentedColormap.from_list(name=r,colors=['white',isau.paperColors(r)]),
                                       bar={'label':label[i], 'ticks':ticks[i], 'orientation':'horizontal'}, interpolation='none', aspect='preserve')
    cb.outline.set_linewidth(0.3)
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    cbs.append(cb)
n = np.array(n)[1:]
ax.set(xlabel='time', xlim=(start,stop), xticklabels=[], ylim=(0,n[-1][1]+11), ylabel='units')
ticks = np.unique(np.concatenate((n.ravel(), n.mean(1))))
ax.set_yticks(ticks, labels=[regs[i//2].upper() if i % 2 else '' for i in range(len(ticks))], rotation=90, va='center')
ticks_labels = ax.get_yticklabels()
ticks = ax.yaxis.get_major_ticks()
for i in range(len(ticks)):
    ticks_labels[i].set_color(isau.paperColors(regs[i//2]) if i % 2 else 'none')
    if i % 2: ticks[i].tick1line.set_visible(False)

# states
height = 0.3; width = 0.8; top = 15 - TOP - height; delta = width + 0.4
for i, state in enumerate(states):
    fma.plotting.plotIntervals(results[f'{state}_interval'].astype(int), ymin=(n[-1][1]+0.5)/(n[-1][1]+11), color=isau.paperColors(state), alpha=0.7)
    cm = Affine2D().scale(1/2.54) + fig.dpi_scale_trans
    rect = Rectangle((LEFT+delta*i,top), width, height, transform=cm, facecolor=isau.paperColors(state), alpha=0.7)
    fig.add_artist(rect)
    fig.text(LEFT+0.1+delta*i, top, state_names[i], fontsize=7, transform=cm, va='bottom')

# color bars
bb = ax.get_position()
width = bb.width / 2
dy = 0.004 # dx = 0.005
y = bb.ymax + dy + 0.001
for cb in reversed(cbs):
    cb.ax.set_position([bb.xmax-width, y, width, dy]) # cb.ax.set_position([x, bb.ymin, dx, bb.height/2])
    cb.ax.set_box_aspect(None)
    cb.ax.tick_params(width=0.3, labelsize=5, pad=0.1, length=1.5)
    y += dy # x += dx


# example zoomed
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig,left,top+height+1, width, 5)
start = results['start_zoom']
stop = results['stop_zoom']
start = 2070
stop = start + 180
n = [[0,0]]
for i, r in enumerate(['nr','hpc','pfc']):
    fr = results[f'fr_zoom_{r}']
    n_units = results[f'units_{r}']
    t = fr[:,0]
    fr = fr[:,1:]
    n.append([n[-1][1], n[-1][1] + n_units])
    im = fma.plotting.plotColorMap(fr.T, x=t, y=n[-1], vmin=0, vmax=.0001, cmap=LinearSegmentedColormap.from_list(name=r,colors=['white',isau.paperColors(r)]),
                                   interpolation='none', aspect='auto')
n = np.array(n)[1:]
ax.set(xlabel='time', xlim=(start,stop), xticks=np.arange(start,stop+1,20), xticklabels=[], ylim=(0,n[-1][1]), ylabel='units')
ticks = np.unique(np.concatenate((n.ravel(), n.mean(1))))
ax.set_yticks(ticks, labels=[regs[i//2].upper() if i % 2 else '' for i in range(len(ticks))], rotation=90, va='center')
ticks_labels = ax.get_yticklabels()
ticks = ax.yaxis.get_major_ticks()
for i in range(len(ticks)):
    ticks_labels[i].set_color(isau.paperColors(regs[i//2]) if i % 2 else 'none')
    if i % 2: ticks[i].tick1line.set_visible(False)


# example zoomed firing rates CHANGE AXES WITH SCALE BARS (NatNeuro Delta induction)
left, top, width, height = isau.pos(ax, fig)
ax = []
for i in range(3):
    ax.append(isau.axes(fig, left+width+1, top+height/3*i, width, height/3))
fr = results[f'fr_zoom_avg']
regs = ('pfc','hpc','nr')
for i, r in enumerate(regs):
    ax[i].plot(fr[:,0], fr[:,i+1], color=isau.paperColors(r))
    ax[i].set(xlim=(start-3,stop), xticks=[], ylim=(0,5), yticks=())
ax[1].set_ylabel('population firing rate (Hz/unit)')
ax[1].yaxis.labelpad = 0.7
ax[-1].yaxis.tick_right(); ax[-1].set(yticks=(0,2.5,5))
ax[-1].spines['right'].set_visible(True); ax[-1].spines["right"].set_position(("outward",5))
ax[-1].set(xlabel='time', xticks=np.arange(start,stop+1,20), xticklabels=[])


fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])