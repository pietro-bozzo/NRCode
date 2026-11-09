import numpy as np
import fmatoolbox as fma
import isautilities as isau
from matplotlib.colors import LinearSegmentedColormap
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig1'
LEFT = 0.7; TOP = 0.5; GAP = 1
axs = {}
regs = ('nr','hpc','pfc')
states = ('sws','rem','other')
state_names = ('nREM','REM','wake')
fig = fma.plotting.figure(n=0, size=(17.5,17.5))


# histology
ax = isau.axes(fig, LEFT, TOP,7.3,2.5)
ax.set(xticks=(), yticks=())
ax.text(.5, .3, 'histology')


# implantation schemes
_ax = isau.axes(fig, GAP, 0,1, 1, ax, 'lwh')
_ax.set(xticks=(), yticks=())
_ax.text(.5, .3, 'headstage')



# example raster
results = isau.load_results(res_root / 'exampleISA.pkl')
ax = isau.axes(fig, 0, GAP,1,5, ax, 'tw')
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
_, top, _, _ = isau.pos(ax)
height = 0.3; width = 0.9; delta = width + 0.2; top -= height + 0.05
for i, state in enumerate(states):
    fma.plotting.plotIntervals(results[f'{state}_interval'].astype(int), ymin=(n[-1][1]+0.5)/(n[-1][1]+11), color=isau.paperColors(state), alpha=0.7)
    isau.rectangle(ax, LEFT+delta*i, top, width, height, facecolor=isau.paperColors(state), alpha=0.7, text=state_names[i])

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
def timescale(ax, start, stop, sec):
    scale = 1 - sec / (stop - start)
    left, top, width, height = isau.pos(ax)
    isau.line(ax, [left+width*scale, left+width], [top+height+0.15, top+height+0.15], color='k', lw=0.6)
    isau.text(ax, left+width*(scale + (1-scale)/2), top+height+0.2, f'{sec} s', fontsize=4, ha='center', va='top')
    return
timescale(ax, start, stop, 500)


# example zoomed
left, top, width, height = isau.pos(ax)
ax = isau.axes(fig, left, top+height+GAP, width, 5)
axs['example_zoom'] = ax
duration = stop - start; old_start = start
start = results['start_zoom']; stop = results['stop_zoom']
# zoom lines
isau.line(ax, [left + width*(start-old_start)/duration, left], [top+height,top+height+1], color='k', lw=0.4, ls='--')
isau.line(ax, [left + width*(stop-old_start)/duration, left + width], [top+height,top+height+1], color='k', lw=0.4, ls='--')
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
ax = []; new_height = height / 6
for i in range(3):
    ax.append(isau.axes(fig, left+width+GAP, top+new_height*i, width, new_height))
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
ax[-1].yaxis.tick_right(); ax[-1].set(yticks=(0,2.5,5))
ax[-1].spines['right'].set_visible(True); ax[-1].spines["right"].set_position(("outward",5))
ax[-1].set(xticks=())
ax[-1].spines.bottom.set_visible(False)

# timescale
timescale(ax[-1], start, stop, 20)


# example cycle
results = isau.load_results(res_root / 'exampleISA.pkl')
ax = isau.axes(fig, LEFT+width+GAP, top+height-2, 3, 2)
axs['cycle_ex'] = ax
duration = stop - start; old_start = start
start = results['start_ex']; stop = results['stop_ex']
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



# FR and CCGs
_, _, _, h = isau.pos(axs['example'])
results = isau.load_results(res_root / 'characterizeRegions.pkl')
fr = results['fr']; ccg_10 = results['ccg_10ms']; ccg_5 = results['ccg_5ms']
width = 2.5; height = 1.2; ax = []; insets = []
for i, r0 in enumerate(ccg_10.reg0.values):
    ax.append([])
    for j, r1 in enumerate(ccg_10.reg0.values):
        ax[-1].append(isau.axes(fig, GAP+i*(width+0.3), j*(height+(h-3*height)/2), width, height, axs['example'], 'l'))
        axs[f'ccg_{r0}_{r1}' if j else f'fr_{r0}'] = ax[-1][-1]
        if j == 1:
            insets.append(isau.axes(fig, width*0.58, -0.1, 0.4, 0.98, ax[-1][-1], 'wh'))
            axs[f'ccg_{r0}_{r1}_inset'] = insets[-1]


conditions = ('nrem','rem','wake'); combinations = (('hpc','nr'),('hpc','pfc'),('nr','pfc'))
for i, r0 in enumerate(ccg_10.reg0.values):
    for j, r1 in enumerate(ccg_10.reg1.values):

        # 1. firing rates
        if j == 0:
            ax[i][j].set_title(f'   {r0.upper()}', size=6, loc='left', pad=-10)
            #b = fma.statistics.hierarchicalBootstrap(tuple(fr[r0].sel(cond=c) for c in conditions), (fr[r0].rat,)*len(conditions))
            fma.plotting.boxPlot(tuple(fr[r0].sel(cond=c) for c in conditions), mode='violin', color=isau.paperColors(conditions), ax=ax[i][j])
            ax[i][j].set(xticks=(), ylim=(0,8), yticks=(0,4,8))

        # 2. CCGs
        else:
            if j == 1:
                r1 = r0; valid = (ccg_10.t.values >= -0.13) & (ccg_10.t.values <= 0.6)
            else:
                r0, r1 = combinations[i]; valid = np.abs(ccg_10.t.values) <= 0.3
            for k, cond in enumerate(conditions):

                ax[i][j].set_title(f'   {r0.upper()}-{r1.upper()}', size=6, loc='left', pad=-10)
                c = ccg_10.sel(cond=cond).sel(reg0=r0,reg1=r1).values
                fma.plotting.semPlot(ccg_10.t[valid], c[:,valid], zscore=2, smooth=1, ax=ax[i][j], color=isau.paperColors(cond))
                if j == 1:
                    c = ccg_5.sel(cond=cond).sel(reg0=r0,reg1=r1).values
                    fma.plotting.semPlot(ccg_5.t, c, zscore=2, smooth=1, ax=insets[i], color=isau.paperColors(cond))
                    if k == 2:
                        ax[i][j].set(xticks=(0,0.2))
                        xl = (-0.06,0.06); yl = np.array(ax[i][j].get_ylim()) * (0.98,0.95)
                        insets[i].set(xlim=xl, xticks=(), ylim=yl, yticks=())
                        insets[i].fill((xl[0],xl[1],xl[1],xl[0]),(yl[0],yl[0],yl[1],yl[1]), fc='white')
                        for s in insets[i].spines:
                            insets[i].spines[s].set_visible(True); insets[i].spines[s].set_linewidth(0.5); insets[i].spines[s].set_color('gray')
                        yl += (0.5,0.4) # compensate for y difference of smaller CCG bin
                        ax[i][j].fill((xl[0],xl[1],xl[1],xl[0]), (yl[0],yl[0],yl[1],yl[1]), fc='none', ec='gray', lw=0.5)
ax[0][0].set(ylabel='firing rate (Hz)')
ax[0][-1].set(xlabel='lag (s)', ylabel='CCG (z-score)')



fig.show()
fma.plotting.saveFigure(fig, output_file, ('png','svg','pdf'))
isau.save_axes_figs(axs, fig, output_file)