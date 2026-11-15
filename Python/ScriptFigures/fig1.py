import numpy as np
import fmatoolbox as fma
import isautilities as isau
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig1'
FW = 17.5; LEFT = 0.7; TOP = 0.5; GAP = 1.5; VSPACE = 1; SAVE = True
axs = {}
regs = ('nr','hpc','pfc')
states = ('sws','rem','other')
state_names = ('nREM','REM','wake')
fig = fma.plotting.figure(n=0, size=(FW,15.5))


# histology
ax = isau.axes(fig, LEFT, TOP,7.3,1.7)
SAVE and ax.imshow(plt.imread(res_root / 'implantation.png'))
bb = ax.get_position()
ax.set_position((0.2/FW, bb.y0, bb.width, bb.height))
ax.axis("off")


# implantation schemes
_ax = isau.axes(fig, -0.5, 0,1, 1, ax, 'lwh')
SAVE and _ax.imshow(plt.imread(res_root / 'HPC-NR-mPFC_allen_new.png'))
_ax.axis("off")



# example raster
results = isau.load_results(res_root / 'exampleISA.pkl')
ax = isau.axes(fig, LEFT, VSPACE,7.3,5, ax, 't', tick_fs=6)
axs['example'] = ax
start = results['start']; stop = results['stop']
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
    _, cb = fma.plotting.plotColorMap(fr.T, x=t, y=n[-1], vmin=0, vmax=2, cmap=LinearSegmentedColormap.from_list(name=r,colors=['white',isau.paperColors(r)]),
                                       bar={'label':label[i], 'ticks':ticks[i], 'orientation':'horizontal'}, interpolation='none', aspect='preserve')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    cbs.append(cb)
n = np.array(n)[1:]
ax.set(xlim=(start,stop), xticks=[], ylim=(0,n[-1][1]+11), ylabel='units')
ax.spines.bottom.set_visible(False)
ticks = np.unique(np.concatenate((n.ravel(), n.mean(1))))
ax.set_yticks(ticks, labels=[regs[i//2].upper() if i % 2 else '' for i in range(len(ticks))], rotation=90, va='center')
ticks_labels = ax.get_yticklabels()
ticks = ax.yaxis.get_major_ticks()
for i in range(len(ticks)):
    ticks_labels[i].set_color(isau.paperColors(regs[i//2]) if i % 2 else 'none')
    if i % 2: ticks[i].tick1line.set_visible(False)

# states
left, top, _, _ = isau.pos(ax)
height = 0.3; width = 0.9; delta = width + 0.2; top -= height + 0.05
for i, state in enumerate(states):
    fma.plotting.plotIntervals(results[f'{state}_interval'].astype(int), ymin=(n[-1][1]+0.5)/(n[-1][1]+11), color=isau.paperColors(state), alpha=0.6)
    isau.rectangle(ax, left+delta*i, top, width, height, facecolor=isau.paperColors(state), alpha=0.6, text=state_names[i])

# color bars
bb = ax.get_position()
width = bb.width / 2
dy = 0.004
y = bb.ymax + dy + 0.001
for cb in reversed(cbs):
    cb.ax.set_position([bb.xmax-width, y, width, dy])
    cb.ax.set_box_aspect(None)
    y += dy

# timescale
def timescale(ax, start, stop, sec, dir='x'):
    left, top, width, height = isau.pos(ax)
    if dir == 'x':
        scale = width * sec / (stop - start)
        x = (left+width-scale, left+width)
        y = (top+height+0.15, top+height+0.15)
        t = (np.mean(x), top+height+0.3)
    else:
        scale = height * sec / (stop - start)
        x = (left+width+0.15, left+width+0.15)
        y = (top+height-scale, top+height)
        t = (left+width+0.3, np.mean(y))
    isau.line(ax, x, y, color='k', lw=0.6)
    isau.text(ax, t[0], t[1], f'{sec} s', fontsize=5, ha='center', va='center')
    return
timescale(ax, start, stop, 500)


# example zoomed
left, top, width, height = isau.pos(ax)
ax = isau.axes(fig, 0, VSPACE, 1, 5, ax, 'tw')
axs['example_zoom'] = ax
duration = stop - start; old_start = start
start = results['start_zoom']; stop = results['stop_zoom']
# zoom lines
isau.line(ax, [left + width*(start-old_start)/duration, left], [top+height,top+height+VSPACE], color='k', ls='--')
isau.line(ax, [left + width*(stop-old_start)/duration, left + width], [top+height,top+height+VSPACE], color='k', ls='--')
n = [[0,0]]
for i, r in enumerate(['nr','hpc','pfc']):
    fr = results[f'fr_zoom_{r}']
    n_units = results[f'units_{r}']
    t = fr[:,0]
    fr = fr[:,1:]
    n.append([n[-1][1], n[-1][1] + n_units])
    fma.plotting.plotColorMap(fr.T, x=t, y=n[-1], vmin=0, vmax=.0001, interpolation='none', aspect='auto',
                                   cmap=LinearSegmentedColormap.from_list(name=r,colors=['white',isau.paperColors(r)]))
n = np.array(n)[1:]
ax.set(xlim=(start,stop), xticks=(), ylim=(0,n[-1][1]), ylabel='units')
ax.spines.bottom.set_visible(False)
ticks = np.unique(np.concatenate((n.ravel(), n.mean(1))))
ax.set_yticks(ticks, labels=[regs[i//2].upper() if i % 2 else '' for i in range(len(ticks))], rotation=90, va='center')
ticks_labels = ax.get_yticklabels()
ticks = ax.yaxis.get_major_ticks()
for i in range(len(ticks)):
    ticks_labels[i].set_color(isau.paperColors(regs[i//2]) if i % 2 else 'none')
    if i % 2: ticks[i].tick1line.set_visible(False)

# timescale
timescale(ax, start, stop, 20)


# example zoomed firing rates
left, top, width, height = isau.pos(ax)
ax = []; new_height = height / 3
for i in range(3):
    ax.append(isau.axes(fig, left+width+GAP, top+new_height*i, width/2, new_height))
    axs[f'example_fr{i}'] = ax[-1]
fr = results['fr_zoom_avg']
regs = ('pfc','hpc','nr')
for i, r in enumerate(regs):
    ax[i].plot(fr[:,0], fr[:,i+1], color=isau.paperColors(r), lw=0.8)
    ax[i].set(xlim=(start-3,stop), xticks=[], ylim=(0,5), yticks=())
    ax[i].spines.left.set_visible(False)
ax[0].set(yticks=[0], yticklabels=()); ax[1].set(yticks=[0], yticklabels=())
ax[1].set_ylabel('population firing rate (Hz/unit)')
ax[1].yaxis.labelpad = 0.7
#ax[-1].yaxis.tick_right(); ax[-1].set(yticks=(0,2.5,5))
#ax[-1].spines['right'].set_visible(True); ax[-1].spines["right"].set_position(("outward",5))
ax[-1].set(xticks=())
ax[-1].spines.bottom.set_visible(False)

# timescale
timescale(ax[-1], start, stop, 20)
timescale(ax[-1], 0, 5, 1, dir='y')


# example cycle
results = isau.load_results(res_root / 'exampleISA.pkl')
_, _, w, h = isau.pos(axs['example_zoom'])
ax = isau.axes(fig, GAP+w*0.6, 0, w*0.4, h*0.33, axs['example_zoom'], 'l')
axs['cycle_ex'] = ax
duration = stop - start; old_start = start
start = results['start_ex']; stop = results['stop_ex']
ax.plot(fr[:,0], fr[:,-1], color=isau.paperColors('nr'), lw=0.8)
ax.set(xlim=(start-3,stop), ylim=(0,5))
ax.axis('off')

ax = isau.axes(fig, GAP+w*0.6, h*0.4, w*0.4, h*0.6, axs['example_zoom'], 'l')
axs['cycle_raster'] = ax
# zoom rectangle
isau.rectangle(ax, LEFT + width*(start-old_start)/duration, top + height*(n[-1,-1]-n[1,0])/n[-1,-1], width*(stop-start)/duration,
               height*n[1,0]/n[-1,-1], fc='none', ec='k', lw=0.4, ls='--')
spikes = results['spike_ex']
fma.plotting.plotIntervals([[start,start+6.8], [stop-1.6,stop]], alpha=0.2)
fma.plotting.plotRaster(spikes, color=isau.paperColors('nr'), compact=True, offset=0, lw=0.8)
ax.set(xlim=(start,stop), xticks=(), ylabel='NR units', ylim=(-.5,results['n_units_nr']-1.5), yticks=())
ax.spines.bottom.set_visible(False); ax.yaxis.labelpad = 0.7
left, top, width, height = isau.pos(ax)
isau.line(fig, [left+width*8/9,left+width], [top+3*height+0.15,top+3*height+0.15], color='k', lw=0.6)
isau.text(fig, left+width*8.5/9, top+3*height+0.2, '2 s', fontsize=4, ha='center', va='top')
perc = 0.35
isau.rectangle(ax, left - width*perc/2, top - height*perc/2, width*(1+perc), height*(1+perc), fc='none', ec='k', lw=0.4, ls='--')
timescale(ax, start, stop, 5)



# FR, BI, and CCGs
_, _, w, h = isau.pos(axs['example']); dh = 2.5; h += dh
results = isau.load_results(res_root / 'characterizeRegions.pkl')
fr = results['fr']; bi = results['bi']; ccg_10 = results['ccg_10ms']; ccg_5 = results['ccg_5ms']
width = w * 0.26; wspace = (w - width*3) / 2; height = h * 0.2; hspace = (h - height*4) / 3; ax = []
for i, r0 in enumerate(ccg_10.reg0.values):
    ax.append([])
    for j in range(4):
        ax[-1].append(isau.axes(fig, GAP+i*(width+wspace), -dh+j*(height+hspace), width, height, axs['example'], 'l'))

conditions = ('nrem','rem','wake'); combinations = (('hpc','nr'),('hpc','pfc'),('nr','pfc'))
for i, r0 in enumerate(ccg_10.reg0.values):
    for j in range(4):

        # 1. firing rates
        if j == 0:
            axs[f'fr_{r0}'] = ax[i][j]
            ax[i][j].set_title(f'   {r0.upper()}', size=6, loc='left', pad=-10)
            # old: non-bootstrapped, ylim: (0,8)
            #fma.plotting.boxPlot(tuple(fr[r0].sel(cond=c) for c in conditions), mode='violin', color=isau.paperColors(conditions), ax=ax[i][j])
            p = fma.statistics.hierarchicalBootstrapTest(tuple(fr[r0].sel(cond=c) for c in conditions), (fr[r0].rat,)*len(conditions), paired=True)
            b = fma.statistics.hierarchicalBootstrap(tuple(fr[r0].sel(cond=c) for c in conditions), (fr[r0].rat,)*len(conditions))
            fma.plotting.boxPlot(b, mode='violin', bw_method=0.5, color=isau.paperColors(conditions), ax=ax[i][j])
            ax[i][j].set(ylim=(1,4))
            fma.plotting.pBar(p[['Condition 1', 'Condition 2', 'Corrected p-value']], draw=(True,False,False,False), ax=ax[i][j])
            ax[i][j].set(xticks=(), ylabel='' if i else 'firing rate (Hz)', ylim=(1,4.5), yticks=(1,4))

        # 2. BI
        elif j == 1:
            axs[f'bi_{r0}'] = ax[i][j]
            ax[i][j].set_title(f'   {r0.upper()}', size=6, loc='left', pad=-10)
            # non-bootstrapped
            fma.plotting.boxPlot(tuple(bi[r0].sel(cond=c) for c in conditions), mode='violin', color=isau.paperColors(conditions), ax=ax[i][j])
            p = fma.statistics.hierarchicalBootstrapTest(tuple(bi[r0].sel(cond=c) for c in conditions), (bi[r0].rat,)*len(conditions), paired=True)
            #b = fma.statistics.hierarchicalBootstrap(tuple(bi[r0].sel(cond=c) for c in conditions), (bi[r0].rat,)*len(conditions))
            #fma.plotting.boxPlot(b, mode='violin', bw_method=None, color=isau.paperColors(conditions), ax=ax[i][j])
            ax[i][j].set(ylim=(-1,1))
            fma.plotting.pBar(p[['Condition 1', 'Condition 2', 'Corrected p-value']], draw=(True,False,False,False), ax=ax[i][j])
            ax[i][j].set(xticks=(), ylabel='' if i else 'burst index', ylim=(-1,1), yticks=(-1,0,1))

        # 3. CCGs
        else:
            ax[i][j].axvline(0, color='k', ls='--', lw=0.5)
            if j == 2:
                r1 = r0; valid = (ccg_5.t.values >= -0.03) & (ccg_5.t.values <= 0.08)
                c = ccg_5.sel(reg0=r0, reg1=r1)
            else:
                r0, r1 = combinations[i]; valid = np.abs(ccg_10.t.values) <= 0.3
                c = ccg_10.sel(reg0=r0, reg1=r1)
            axs[f'ccg_{r0}_{r1}'] = ax[i][j]
            ax[i][j].set_title(f'   {r0.upper()}-{r1.upper()}', size=6, loc='left', pad=-10)
            for k, cond in enumerate(conditions):
                fma.plotting.semPlot(c.t[valid], c.sel(cond=cond)[:,valid], zscore=2, smooth=1, ax=ax[i][j], color=isau.paperColors(cond))
            ax[i][j].set(ylabel='' if i else 'CCG (z-score)')
[ax[i][-1].set(xlabel='lag (s)') for i in range(3)]



fig.show()
SAVE and fma.plotting.saveFigure(fig, output_file, ('png','svg','pdf'))
SAVE and isau.save_axes_figs(axs, fig, output_file)