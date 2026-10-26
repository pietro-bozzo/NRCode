import fmatoolbox as fma
import isautilities as isau
import pathlib
import numpy as np

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig5'
LEFT = 0.8; TOP = 0.5; VSPACE = 1.2
fig = fma.plotting.figure(n=0, size=(17.5,10.9))
axs = {}

# sleep2 vs sleep 1 results
results = isau.load_results(res_root / 'ISACouplHpcPfc.pkl')
events = results['events']; rate = results['rate']
width = 2
ax = [isau.axes(fig, LEFT+(0.5+width)*i, TOP, width, 3) for i in range(3)]
for n, a in zip(('rate_incr_rip','rate_incr_delta','rate_incr_spin'), ax): axs[n] = a

# input zeros with animal means (computed excluding zeros)
def fill_group(x):
    x_masked = x.where(x != 0) # zeros -> NaN, per rat
    group_mean = x_masked.mean(dim='rat', skipna=True) # mean excluding zeros, per rat
    return x.where(x != 0, group_mean)                 # broadcast fill back over this group's rat entries
rate_inputted = rate.groupby('rat').map(fill_group)


# rate increase
ylim = ((-0.5,1.5), (-0.5,1.5), (-0.5,1.5))
yticks = ((-0.5,0,0.5,1,1.5),) * 3
for i, name in enumerate(events):
    ax[i].axhline(0, color='#aaaaaa', ls='--', lw=0.6)
    data = np.log(rate_inputted.sel(coupl=False,event=name,epoch='sleep1') / rate_inputted.sel(coupl=False,event=name,epoch='sleep0'))
    isau.pair_violin(tuple(data.sel(isa=a) for a in [True,False]), (data.rat,)*2, color=isau.paperColors((name,'shuffle')),
                      mode='zero', label=('ISA','nISA'), ax=ax[i], ylim=ylim[i], yticks=yticks[i])
ax[0].set(ylabel='rate increase (log)')


# couples' rate increase
left, top, width, height = isau.pos(ax[-1], fig)
ax = [isau.axes(fig, left+width+2+(0.5+width)*i, top, width, 3) for i in range(3)]
for n, a in zip(('coupl_rate_incr_rip','coupl_rate_incr_delta','coupl_rate_incr_spin'), ax): axs[n] = a
ylim = ((-0.5,1.5), (-0.5,1.5), (-0.5,1.5))
yticks = ((-0.5,0,0.5,1,1.5),) * 3
for i, name in enumerate(events):
    ax[i].axhline(0, color='#aaaaaa', ls='--', lw=0.6)
    data = np.log(rate_inputted.sel(coupl=True,event=name,epoch='sleep1') / rate_inputted.sel(coupl=True,event=name,epoch='sleep0'))
    isau.pair_violin(tuple(data.sel(isa=a) for a in [True,False]), (data.rat,)*2, color=isau.paperColors((name,'shuffle')),
                      mode='zero', label=('ISA','nISA'), ax=ax[i], ylim=ylim[i], yticks=yticks[i])
ax[0].set_ylabel('rate increase of\ncoupled events (log)')




fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])
isau.save_axes_figs(axs, fig, output_file)