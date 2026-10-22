import fmatoolbox as fma
import isautilities as isau
import pathlib
import matplotlib.ticker

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig3'
LEFT = 0.8; TOP = 0.3; VSPACE = 0.8
axs = {}
fig = fma.plotting.figure(n=0, size=(17.5,10.9))

# example
ax = isau.axes(fig, LEFT, TOP, 7, 3)
axs['example'] = ax

# oscillations are more likely in ISA
results = isau.load_results(res_root / 'ISACouplHpcPfc.pkl')
left, top, width, height = isau.pos(ax, fig)
ax = [isau.axes(fig, left+width+1.5+i*(0.4+2.6), top, 2.4, height) for i in range(3)]
for n, a in zip(('base_rate_rip','base_rate_delta','base_rate_spin'), ax): axs[n] = a
events = results['events']; rate = results['rate']; perc = results['perc']
ylim = ((0,0.6), (0,0.6), (0,0.6))
yticks = ((0,0.2,0.4,0.6), (0,0.2,0.4,0.6), (0,0.2,0.4,0.6))
for i, e in enumerate(events):
    tc = rate.sel(event=e, epoch='sleep0', coupl=False)
    isau.pair_violin(tuple(tc.sel(isa=a) for a in [True,False]), (tc.rat,)*2, color=isau.paperColors((e,'shuffle')), label=('ISA','nISA'),
                     dx=0.3, ax=ax[i], ylim=ylim[i], yticks=yticks[i],)
ax[0].set(ylabel='rate of events (Hz)')


# ISA entrains oscillations
results = isau.load_results(res_root / 'HpcPfcEntrainment.pkl')
_, top, width, height = isau.pos(ax[0], fig)
new_width = 3
ax = [isau.axes(fig, LEFT+i*(0.1+new_width), top+height+VSPACE, new_width, 2) for i in range(3)]
for n, a in zip(('entrain_rip','entrain_delta','entrain_spin'), ax): axs[n] = a
tc = results['tc']; kur = results['kur']
ylim = ((0.25,1.75), (0.25,1.75), (0.25,1.75)) # ((0.5,1.5), (0.45,1.55), (0.25,1.75))
yticks = ((0.25,1,1.75), (), ())
for i, e in enumerate(events):
    fma.plotting.semPlot(tc.phi, tc.sel(ev=e,shuf=True), saxis=2, polar=2, color='k', ax=ax[i], smooth=1, aprop={'facecolor':'none','lw':0.6})
    fma.plotting.semPlot(tc.phi, tc.sel(ev=e,shuf=False), polar=2, color=isau.paperColors(e), ax=ax[i], smooth=1.5)
    ax[i].set(ylim=ylim[i], yticks=yticks[i])
    ax[i].spines['left'].set_visible(i==0)
isau.ISAxes(ax, double=True)
ax[0].set(ylabel='tuning curve')


# R
left, top, width, height = isau.pos(ax[-1], fig)
new_width = 1.3
ax = [isau.axes(fig, left+width+1.5+i*(0.7+new_width), top, new_width, height) for i in range(3)]
for n, a in zip(('R_rip','R_delta','R_spin'), ax): axs[n] = a
ylim = ((0,0.25), (0,0.25), (0,0.25))
yticks = ((0,0.1,0.2,0.3), (0,0.1,0.2,0.3), (0,0.1,0.2,0.3))
for i, e in enumerate(events):
    tc = kur.sel(ev=e, val='R')
    isau.pair_violin((tc.sel(shuf=False), tc.isel(shuf=17)), (tc.rat,)*2, color=(isau.paperColors((e)),'k'), label=('obs','sur'),
                     vprop=({},{'facecolor':'none','lw':0.6}), ax=ax[i], ylim=ylim[i], yticks=yticks[i])
ax[0].set(ylabel='modulation strength R')



# coupled oscillations are more likely in ISA
left, top, width, height = isau.pos(ax[-1], fig)
new_width = 1.3
ax = [isau.axes(fig, LEFT+i*(0.7+new_width), top+height+VSPACE+0.2, new_width, height) for i in range(3)]
for n, a in zip(('coupl_rate_rip','coupl_rate_delta','coupl_rate_spin'), ax): axs[n] = a
ylim = ((0,0.035), (0,0.015), (0,0.08))
yticks = ((0,0.01,0.02,0.03), (0,0.01,0.02), (0,0.05,0.1))
for i, e in enumerate(events):
    tc = rate.sel(event=e, epoch='sleep0', coupl=True)
    isau.pair_violin(tuple(tc.sel(isa=a) for a in [True,False]), (tc.rat,)*2, color=isau.paperColors((e,'shuffle')), label=('ISA','nISA'),
                     ax=ax[i], ylim=ylim[i], yticks=yticks[i])
    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1,-1))
    ax[i].yaxis.set_major_formatter(formatter)
    ax[i].yaxis.get_offset_text().set_fontsize(4)
ax[0].set(ylabel='rate of coupled events (Hz)')


# P(coupling) EXCLUDED
# left, top, width, height = isau.pos(ax[-1], fig)
# new_width = 1.3
# ax = [isau.axes(fig, left+width+1+i*(0.7+new_width), top, new_width, height) for i in range(3)]
# ylim = ((0,0.1), (0,0.5), (0,0.5))
# yticks = ((0,0.01,0.02,0.03), (0,0.01,0.02), (0,0.05,0.1))
# for i, e in enumerate(events):
#     tc = perc.sel(event=e, epoch='sleep0', binomial='p')
#     isau.pair_violin(tuple(tc.sel(isa=a) for a in [True,False]), (tc.rat,)*2, color=isau.paperColors((e,'shuffle')), label=('ISA','nISA'),
#                      ax=ax[i], ylim=ylim[i], yticks=yticks[i])
# ax[0].set(ylabel='P of coupling')



fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])
isau.save_axes_figs(axs, fig, output_file)