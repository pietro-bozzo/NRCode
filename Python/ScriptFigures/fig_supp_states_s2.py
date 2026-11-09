import numpy as np
import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'FigSuppStatesS2'
LEFT = 0.9; TOP = 0.5
axs = {}
fig = fma.plotting.figure(n=0, size=(14.5,4.5))
results = isau.load_results(res_root / 'ISAStates.pkl')


# % of states spent in ISA
ax = isau.axes(fig, LEFT, TOP, 2.1, 2.3)
axs['states_inv'] = ax
states = results['states']; f = results['fraction'].sel(val='reverse'); f2 = results['fraction_s2'].sel(val='reverse')
isau.pair_violin((f.sel(state='sws')*100, f2.sel(state='sws')*100), (f.rat,f2.rat), ax=ax, ylim=(0,80), pdraw=(True,)*4,
                 color=isau.paperColors(('sws','sws')), label=('baseline','post-task'), ylabel='nREM time spent in ISA (%)')


# % of ISA spent in states
left, top, width, height = isau.pos(ax)
rect_space = 0.5; rect_dim = 0.3
ax = isau.axes(fig, left+width+1.2, top+rect_space, height-rect_space, height-rect_space)
left, _, width, _ = isau.pos(ax)
axs['states'] = ax
medians = []
for i, state in enumerate(states):
    fd = results['fraction_s2'].sel(state=state)
    medians.append([np.nanmedian(fd.sel(val='ref') * 100), np.nanmedian(fd.sel(val='isa') * 100)])
    ax.scatter(fd.sel(val='ref')*100, fd.sel(val='isa')*100, s=20, alpha=0.5, color=isau.paperColors(state), edgecolor='none')
ax.plot(ax.get_xlim(), ax.get_ylim(), color='gray', ls='--', lw=0.5, zorder=-1)
ax.set(xlabel='time spent in state (%)', xlim=(ax.get_xlim()[0],80), xticks=(0,40,80),
       ylabel='ISA time spent in state (%)', ylim=(ax.get_ylim()[0],105), yticks=(0,50,100))
widths = [0,0]
for i, state in enumerate(states):
    m1, m2 = medians[i]
    isau.rectangle(ax, left+widths[0], top, width*m1/100, rect_dim, ec=isau.paperColors(state),
                   fc=isau.paperColors(state,'70'), lw=0.5, text=f'{m1:.0f}')
    isau.rectangle(ax, left+width+rect_space-rect_dim, top+rect_space+widths[1], rect_dim, width*m2/100, ec=isau.paperColors(state),
                   fc=isau.paperColors(state,'70'), lw=0.5, text=f'{m2:.0f}')
    widths[0] += width*m1/100; widths[1] += width*m2/100

# nREM interval duration
left, top, width, height = isau.pos(ax)
ax = isau.axes(fig, left+width+1.7, top-rect_space, 2.5, height+rect_space)
axs['nrem_dur'] = ax
f = results['duration_s2']['sws']
fma.plotting.semPlot(f.dur/60, f*100, color=isau.paperColors('sws'), ci='nansem')
stat = fma.statistics.clusterPermutationTest(f*100,50)
ax.set(xlabel='state interval duration (min)',ylim=[0,75])
fma.plotting.pHorzLine(stat['sig_mask'], f.dur/60, color='k')
ax.set(xlabel='nREM interval duration (min)', ylabel='time spent in ISA (%)', xlim=(ax.get_xlim()[0],10))


# nREM interval dynamics
f = results['interval_s2'].sel(state='sws')
left, top, width, height = isau.pos(ax)
rect_space = 0.5; rect_dim = 0.3
ax = isau.axes(fig, left+width+0.85, top+rect_space, 1, height-rect_space)
axs['nrem_elapsed0'] = ax
fma.plotting.semPlot(f.win, f.sel(val='start')*100, color=isau.paperColors('sws'))
ax.axvline(50, ls='--', lw=0.6, color='k')
ax.set(xlabel='                  elapsed nREM time (s)', xlim=(ax.get_xlim()[0],50), xticks=(0,50),
       xticklabels=(0,'50    '), ylim=(0,60))

left, top, width, height = isau.pos(ax)
ax = isau.axes(fig, left+width+0.2, top, 1, height)
axs['nrem_elapsed1'] = ax
fma.plotting.semPlot(f.win,f.sel(val='stop')*100,color=isau.paperColors('sws'))
ax.axvline(0, ls='--', lw=0.6, color='k')
ax.set(xlim=(0, ax.get_xlim()[1]), xticks=(0,50), xticklabels=('    -50',0), ylim=(0,60), yticks=())
ax.spines['left'].set_visible(False)

cumw = 0
for w, c, t in zip((0.2,0.6,0.2), isau.paperColors(('shuffle','nr','shuffle')), (None,'ISA',None)):
    isau.rectangle(ax, left+cumw, top-rect_space, (2*width+0.2)*w, rect_dim, ec=c, fc=c[:-2]+'70', lw=0.5, text=t)
    cumw += (2*width+0.2)*w



fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])
isau.save_axes_figs(axs, fig, output_file)