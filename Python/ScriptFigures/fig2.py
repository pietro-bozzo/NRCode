import numpy as np
import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig2'
LEFT = 0.7; TOP = 0.3
axs = {}
fig = fma.plotting.figure(n=0, size=(17.5,10.9))


# ISA is exclusive to NR
results = isau.load_results(res_root / 'detectISA.pkl')
ax = isau.axes(fig, LEFT, TOP, 4, 3)
axs['unique'] = ax
labels = ['NR','TH','iHPC','vHPC','dHPC','PFC','V1','AMY']
regs = results['regs']
color = isau.paperColors([regs[i] for i in [0,0,2,2,2,5,5,7]], ['ff','90'] + ['ff']*(len(regs)-2))
b = results['ISA_fraction']
p = results['p_val']
fma.plotting.boxPlot(b, mode='violin', color=color, widths=1, label=labels)
fma.plotting.pBar(p[["Condition 2","Condition 2","Corrected p-value"]], x=range(len(regs)))
ax.set(ylabel='ISA time (%)', ylim=(0,40))
ax.tick_params(axis='x', pad=4, rotation=45)


# example cycle
results = isau.load_results(res_root / 'exampleISA.pkl')
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+0.4, top, 2.6, height)
axs['cycle_ex'] = ax
start = results['start_ex']; stop = results['stop_ex']
spikes = results['spike_ex']
fma.plotting.plotIntervals([[start,start+6.8], [stop-1.6,stop]], alpha=0.2)
fma.plotting.plotRaster(spikes, color=isau.paperColors('nr'), compact=True, offset=0, lw=0.8)
ax.set(xlabel='time', xlim=(start,stop), xticks=range(start,stop+1,2), xticklabels=(), ylabel='NR units', ylim=(-.5,results['n_units_nr']-1.5), yticks=())
ax.yaxis.labelpad = 0.7


# typical cycle
results = isau.load_results(res_root / 'regionStatusFiring.pkl')
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+1.2, top, 3.5, height)
axs['cycle_avg'] = ax
fr_on = results['fr_on']; fr_off = results['fr_off']
for i, reg in enumerate(['nr']):
    x = np.concatenate((fr_on.phi,fr_off.phi))
    data = np.concatenate((fr_on.sel(reg=reg),fr_off.sel(reg=reg)),axis=1)
    fma.plotting.semPlot(x,data,polar=2,smooth=1.5,ax=ax,color=isau.paperColors(reg))
ax.set(ylabel='population firing\n(Hz/neurons)', ylim=(0.5,2.5), yticks=(0.5,1.5,2.5))
isau.ISAxes(ax, double=True)


# durations
results = isau.load_results(res_root / 'characterizeISA.pkl')
isa, on, off, cycles = results['isa'], results['on'], results['off'], results['cycles']

left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+1.2, top, 3, height*0.4)
axs['cycle_pdf'] = ax
y = []
for i, ra in enumerate(cycles):
    x, _y = fma.analysis.PDF(np.diff(cycles[ra]['nr']), bins=np.linspace(0,30,40))
    y.append(_y)
fma.plotting.semPlot(x,y,color=isau.paperColors('nr'))
ax.set(xlim=(0,30),xlabel='ISA cycle duration (s)',ylim=(ax.get_ylim()[0],0.1))#,yticks=np.linspace(0,0.1,4))

left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left, TOP+3-height, 3, height)
axs['epoch_pdf'] = ax
y = []
for i, ra in enumerate(isa):
    x, _y = fma.analysis.PDF(np.diff(isa[ra]['nr']), bins=np.linspace(0,300,40))
    y.append(_y)
fma.plotting.semPlot(x,y,color=isau.paperColors('nr'))
ax.set(xlim=(0,300),xlabel='ISA epoch duration (s)',ylim=(ax.get_ylim()[0],0.0085))#,yticks=np.linspace(0,0.1,4))
ax.set_ylabel('                       normalized count',labelpad=3)




# % of states spent in ISA
results = isau.load_results(res_root / 'ISAStates.pkl')
ax = isau.axes(fig, LEFT+0.2, TOP+4, 2.1, 2.3)
axs['states_inv'] = ax
states, f = results['states'], results['fraction'].sel(val='reverse')
b = fma.statistics.hierarchicalBootstrap(tuple(f.sel(state=s)*100 for s in states), tuple(f.rat for s in states))
p = fma.statistics.hierarchicalBootstrapTest(tuple(f.sel(state=s) for s in states), tuple(f.rat for s in states), paired=True)
fma.plotting.boxPlot(b, mode='violin', color=isau.paperColors(states), label=('nREM','REM','wake'))
ax.set(ylim=(ax.get_ylim()[0],65))
fma.plotting.pBar(p[['Condition 1', 'Condition 2', 'Corrected p-value']])
ax.set(ylabel='state time spent in ISA (%)', xlim=(-0.5,2.5), ylim=(ax.get_ylim()[0],80))
ax.tick_params(axis='x', pad=10)


# % of ISA spent in states
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+1.2, top, height, height)
axs['states'] = ax
for i, state in enumerate(states):
    fd = results['fraction'].sel(state=state)
    ax.scatter(fd.sel(val='ref')*100, fd.sel(val='isa')*100, s=20, alpha=0.5, color=isau.paperColors(state))
ax.plot(ax.get_xlim(), ax.get_ylim(), color='gray', ls='--', label='1:1', zorder=-1)
ax.set(xlabel='time spent in state (%)', xlim=(ax.get_xlim()[0],80), xticks=(0,40,80),
       ylabel='ISA time spent in state (%)', ylim=(ax.get_ylim()[0],105), yticks=(0,50,100))


# nREM interval duration
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+1.3, top, 2.5, height)
axs['nrem_dur'] = ax
f = results['duration']['sws']
fma.plotting.semPlot(f.dur/60, f*100, color=isau.paperColors('sws'), ci='nansem')
stat = fma.analysis.clusterPermutationTest(f*100,50)
ax.set(xlabel='state interval duration (min)',ylim=[0,60])
fma.plotting.pHorzLine(stat['sig_mask'], f.dur/60, color='k')
ax.set(xlabel='nREM interval duration (min)', ylabel='time spent in ISA (%)', xlim=(ax.get_xlim()[0],10))


# nREM interval dynamics
f = results['interval'].sel(state='sws')
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+0.85, top, 1, height)
axs['nrem_elapsed0'] = ax
fma.plotting.semPlot(f.win, f.sel(val='start')*100, color=isau.paperColors('sws'))
ax.axvline(50, ls='--', lw=0.6, color='k')
ax.set(xlabel='                  elapsed nREM time (s)', xlim=(ax.get_xlim()[0],50), ylim=(0,60))

left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+0.2, top, 1, height)
axs['nrem_elapsed1'] = ax
fma.plotting.semPlot(f.win,f.sel(val='stop')*100,color=isau.paperColors('sws'))
ax.axvline(0, ls='--', lw=0.6, color='k')
ax.set(xlim=(0, ax.get_xlim()[1]), xticks=(0,50), xticklabels=(-50,0), ylim=(0,60), yticks=())
ax.spines['left'].set_visible(False)


# concatenated nREM
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+0.85, top, 2.5, height)
axs['nrem_concat'] = ax
f = results['shifted'].sel(state='sws')
fma.plotting.semPlot(f.win/60, f*100, color=isau.paperColors('sws'), ci='nansem')
stat = fma.analysis.clusterPermutationTest(f*100,57)
ax.set(ylim=(0,60))
fma.plotting.pHorzLine(stat['sig_mask'], f.win/60, color='k', dy=0)
ax.set(xlabel='concatenated nREM time (min)', xlim=(ax.get_xlim()[0],80))




# MAYBE: thinner line On OFF, a bit more ticks padding
results = isau.load_results(res_root / 'unitEntrainment.pkl')
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, LEFT, top+height+1, 2.5, 2.5)
axs['unit_colormap'] = ax
tc = results['tc'].sel(reg='nr')
kur = results['kur'].sel(reg='nr',val='psi')

# separate in (non-)significantly entrained and sort separately
signif = kur.h.values
n_units = [sum(~signif), sum(signif)]
data_signif = tc.isel(reg=signif).values
order = np.argsort(kur.isel(reg=signif)).values
data_signif = data_signif[order]
data_nsignif = tc.isel(reg=~signif).values
order = np.argsort(kur.isel(reg=~signif)).values
data_nsignif = data_nsignif[order]
data = np.concatenate((np.full((200,tc.shape[1]),np.nan), data_nsignif, np.full((200,tc.shape[1]),np.nan), data_signif))

fma.plotting.plotColorMap(data, aspect='auto', x=tc.phi, vmin=-1, vmax=1, zscore='all')
ax.plot([6.6]*2, [200,200+n_units[0]], isau.paperColors('shuffle'), [6.6]*2, [400+n_units[0],400+sum(n_units)], isau.paperColors('nr'), lw=2)
ax.spines['left'].set_visible(False)
ax.set_yticks([200+n_units[0]/2,400+n_units[0]+n_units[1]/2], labels=n_units, rotation=90, va='center')
ax.set(ylabel='units'); isau.ISAxes(ax)


# preferred phase
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+1.3, top, 2.8, height)
axs['status'] = ax
psi = np.mod(kur, 2*np.pi)
y = []
for rat in np.unique(psi.rat):
    x, _y = fma.analysis.PDF(psi.sel(rat=rat).isel(reg=psi.h.sel(rat=rat)), mode='polar', bandwidth=0.03)
    y.append(_y)
fma.plotting.semPlot(x, y, color=isau.paperColors('nr'))
ax.set(ylabel='preferred status\ndistribution')
isau.ISAxes(ax, double=True)


# prefer ON / OFF proportions
left, top, width, height = isau.pos(ax, fig)
ax = isau.axes(fig, left+width+1.1, top, 0.8, height)
axs['on_pref'] = ax
import statsmodels.stats.proportion as smsp
limits = [[-1,10], [np.pi/4,3*np.pi/4], [5*np.pi/4,7*np.pi/4]]
alpha = ('e0','70')
count = []; nobs = []; color = []
for k in range(2):
    this_data = psi.h.values & (psi.values > limits[k][0]) & (psi.values < limits[k][1])
    count.append(this_data.sum())
    nobs.append(len(this_data))
    color.append(isau.paperColors('nr', alpha[k]))
low, high = smsp.proportion_confint(count=count, nobs=nobs, alpha=0.05, method='wilson')
p = np.array(count) / np.array(nobs)
low = p - low
low[low<0] = 0
high = high - p
ax.bar(range(2), np.array(p)*100, yerr=np.stack((low,high),axis=1)*100, color=color)
ax.set(xlim=(-0.8,1.8), ylabel='% of entrained units', ylim=(0,100))
ax.set_xticks(range(2), labels=('all','ON'), fontsize=7)

# avalanche cartoon ? Not again
results = isau.load_results(res_root / 'criticalExponents.pkl')
left, top, width, height = isau.pos(ax, fig)
new_height = 0.85 * height / 2
ax = [isau.axes(fig, left+width+1.1, top, 2.8, new_height)]
for n, a in zip(('aval_S','aval_D'), ax): axs[n] = a
ax.append(isau.axes(fig, left+width+1.1, top+height-new_height, 2.8, new_height))

# S, D
# FIND EXAMPLE WITH BEST range diff: 2? 7? 18? 22?
fit = results['fit_states']
pdf = results['pdf_states']
def plotPowerLaw(pdf, fit, ax, name=None, label=None, **kwargs):
    fma.plotting.plotXY(pdf, label=label, ax=ax, lw=1, **kwargs)
    ax.axvline(fit.sel(val='xmin'), ls='--', lw=0.8, **kwargs)
    ax.axvline(fit.sel(val='xmax'), ls='--', lw=0.8, **kwargs)
    ax.set(xscale='log', yscale='log')
    if name is not None:
        ax.set(xlabel=name, ylabel=f'p({name})')
    ax.set_xticks((1,10,100))
    return
session = 2; i = 0; s = list(pdf[session].keys())[i] # nREM
regs = []
for j, r in enumerate(pdf[session][s]):
    regs.append(r)
    color = isau.paperColors(r)
    plotPowerLaw(pdf[session][s][r]['dur'], fit.sel(reg=r,state=s,q='D')[session], name='D', color=color, ax=ax[1])
    plotPowerLaw(pdf[session][s][r]['size'], fit.sel(reg=r,state=s,q='S')[session], name='S', color=color, ax=ax[0])

# criticality in nREM
exponents = results['exp_states']
left, top, width, _ = isau.pos(ax[0], fig)
ax = isau.axes(fig, left+width+1.1, top, width, height)
axs['dcc'] = ax
dcc = []; #label = []
for i, r in enumerate(regs):
    valid = np.all(~np.isnan(fit.sel(reg=r,state='sws',val='p')), axis=1) # sessions where a fit was performed
    h = fit.sel(reg=r,state='sws',val='p')[valid].values < 0.05
    signif = np.all(h, axis=1) # sessions where the fit was successful
    #label.append(f'{signif.mean():.2f}') # percentage of successful fits
    dcc.append(exponents.sel(reg=r,state='sws',val='DCC'))
    dcc[-1][~valid] = np.nan
    dcc[-1][valid][~signif] = np.nan
valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) | ~np.isnan(dcc[2]) # sessions where at least one condition was significant
for j in range(len(dcc)):
    dcc[j] = dcc[j][valid]
    dcc[j][np.isnan(dcc[j])] = dcc[j][~np.isnan(dcc[j])].max() # input non significant fits with highest observed distance
b = fma.statistics.hierarchicalBootstrap(tuple(dcc), tuple(d.rat for d in dcc))
p = fma.statistics.hierarchicalBootstrapTest(tuple(dcc), tuple(d.rat for d in dcc), paired=True)
fma.plotting.boxPlot(b, mode='violin', color=isau.paperColors(regs), label=('NR','HPC','PFC'))
ax.set(ylim=(0,0.5))
fma.plotting.pBar(p[['Condition 1', 'Condition 2', 'Corrected p-value']])
ax.set(xlim=(-0.5,2.5), ylabel='DCC', ylim=(0,0.55))



fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])
isau.save_axes_figs(axs, fig, output_file)