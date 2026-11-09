import fmatoolbox as fma
import isautilities as isau
import numpy as np
import xarray as xr
import pandas as pd
import statsmodels.formula.api as smf
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig5'
FW = 17.5; FH = 10.9; LEFT = 1; TOP = 0.5; GAP = 1
fig = fma.plotting.figure(n=0, size=(FW,FH))
axs = {}


# ISA proportion is unchanged
results = isau.load_results(res_root / 'detectISA.pkl')
ax = isau.axes(fig, FW/2+GAP/2, TOP, 1.95, 2.5)
axs['ISA_frac_s2'] = ax
frac1, frac2 = results['ISA_s1_s2']; frac1 = frac1.sel(region='nr')
isau.pair_violin((frac1*100,frac2*100), (frac1.rat,frac2.rat), ax=ax, color=isau.paperColors(('nr','nr')),
                 label=('baseline','post-task'), ylabel='ISA time (%)', ylim=(0,37), pdraw=(True,)*4)
ax.set(ylim=(0,40))

# ISA period is unchanged
results = isau.load_results(res_root / 'characterizeISA.pkl')
cycles1 = results['cycles']; cycles2 = results['cycles_s2']
# sleep 1
left, top, width, height = isau.pos(ax)
ax = isau.axes(fig, left+width+GAP, top, width, height*0.4)
axs['ISA_cycle_s1'] = ax
data = [[], []]; rat = [[], []]; y = []
for i, ra in enumerate(cycles1):
    data[0].append(np.diff(cycles1[ra]['nr']).ravel())
    rat[0].append(np.full_like(data[0][-1], ra, dtype=int))
    x, _y = fma.analysis.PDF(data[0][-1], bins=np.linspace(0,30,40))
    y.append(_y)
fma.plotting.semPlot(x, y, color=isau.paperColors('nr'))
ax.set(xlim=(0,30),ylim=(ax.get_ylim()[0],0.1))
# sleep 2
ax = isau.axes(fig, left+width+GAP, top+height*0.6, width, height*0.4)
axs['ISA_cycle_s2'] = ax
y = []
for i, ra in enumerate(cycles2):
    data[1].append(np.diff(cycles2[ra]['nr']).ravel())
    rat[1].append(np.full_like(data[1][-1], ra, dtype=int))
    x, _y = fma.analysis.PDF(data[1][-1], bins=np.linspace(0,30,40))
    y.append(_y)
fma.plotting.semPlot(x,y,color=isau.paperColors('nr'))
ax.set(xlim=(0,30),xlabel='ISA cycle duration (s)',ylim=(ax.get_ylim()[0],0.1))
ax.set_ylabel('                  normalized count\n                   baseline       post-task',labelpad=2)
# stats
p = fma.statistics.hierarchicalBootstrapTest((np.concatenate(data[0]), np.concatenate(data[1])),
                                             (np.concatenate(rat[0]), np.concatenate(rat[1])), paired=True)
line_x = left + width + GAP + width * 0.85
isau.line(ax, (line_x,line_x), (top+0.25,top+1.8), color='k', lw=0.5)
isau.text(ax, left+width+line_x+0.05, top+0.8, 'n. s.', fontsize=5)


# nREM at beginning of sleep session
results = isau.load_results(res_root / 'ISAStates.pkl')
left, _, width, _ = isau.pos(ax)
ax = isau.axes(fig, left+width+GAP, top, width, height)
axs['nrem_concat'] = ax
f = results['shifted'].sel(state='sws'); f2 = results['shifted_s2'].sel(state='sws')
fma.plotting.semPlot(f.win/60, f*100, color=isau.paperColors('sws'), ci='nansem')
fma.plotting.semPlot(f2.win/60, f2*100, color=isau.paperColors('sws','50'), ci='nansem')
stat = fma.statistics.clusterPermutationTest(f,f2)
ax.set(ylim=(0,70))
fma.plotting.pHorzLine(stat['sig_mask'], f.win/60, color='k', dy=0)
ax.set(xlabel='concatenated nREM\ntime (min)', xlim=(ax.get_xlim()[0],80), ylabel='time spent in ISA (%)')



# sleep2 vs sleep 1 results
results = isau.load_results(res_root / 'ISACouplHpcPfc.pkl')
events = results['events']; rate = results['rate']
width = 2
left, top, _, height = isau.pos(ax)
ax = [isau.axes(fig, LEFT+(0.5+width)*i, top+height+GAP, width, 3) for i in range(3)]
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
left, top, width, height = isau.pos(ax[-1])
ax = [isau.axes(fig, FW/2+GAP/2+(0.5+width)*i, top, width, 3) for i in range(3)]
for n, a in zip(('coupl_rate_incr_rip','coupl_rate_incr_delta','coupl_rate_incr_spin'), ax): axs[n] = a
ylim = ((-0.5,1.5), (-0.5,1.5), (-0.5,1.5))
yticks = ((-0.5,0,0.5,1,1.5),) * 3
for i, name in enumerate(events):
    ax[i].axhline(0, color='#aaaaaa', ls='--', lw=0.6)
    data = np.log(rate_inputted.sel(coupl=True,event=name,epoch='sleep1') / rate_inputted.sel(coupl=True,event=name,epoch='sleep0'))
    isau.pair_violin(tuple(data.sel(isa=a) for a in [True,False]), (data.rat,)*2, color=isau.paperColors((name,'shuffle')),
                      mode='zero', label=('ISA','nISA'), ax=ax[i], ylim=ylim[i], yticks=yticks[i])
ax[0].set_ylabel('rate increase of\ncoupled events (log)')



# PFC criticality
_, top, _, height = isau.pos(ax[-1])
ax = isau.axes(fig, LEFT, top+height, 2.5, height)
axs['dcc_pfc_s2'] = ax
results = isau.load_results(res_root / 'criticalExponents.pkl')
fit1 = results['fit_ISA']; exp1 = results['exp_ISA']; fit2 = results['fit_ISA_s2']; exp2 = results['exp_ISA_s2']
base = fit2.base; states = fit2.state.values; regs = fit2.reg.values
fit = xr.concat((fit1,fit2), dim='sleep', join='outer'); fit = fit.sel(base=base)
exp = xr.concat((exp1,exp2), dim='sleep', join='outer'); exp = exp.sel(base=base)
for i, r in enumerate(regs):
    dcc = []; valid = np.array([False], dtype=bool)
    for sl in range(2):
        for k, s in enumerate(states):
            h = fit.sel(sleep=sl, reg=r, state=s, val='p').values < 0.05 # failed sessions get h False since p is nan
            signif = np.all(h, axis=1) # sessions where the fit was performed and significant for S and D
            dcc.append(exp.sel(sleep=sl, reg=r, state=s, val='DCC').copy(deep=True))
            dcc[-1][~signif] = np.nan
            valid = np.logical_or(valid, signif) # sessions where at least one condition was significant

        #valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) # sessions where at least one condition was significant
    rat = []; session = []; isa = []; time = []
    relable = lambda x : np.unique(x, return_inverse=True)[1]
    for k in range(len(dcc)):
        dcc[k] = dcc[k][valid]
        dcc[k][np.isnan(dcc[k])] = dcc[k][~np.isnan(dcc[k])].max() # input non significant fits with highest observed distance
        rat.append(dcc[k].rat)
        session.append(relable(dcc[k].base)) # relabel sessions in {0, N}
        isa.append(dcc[k].state)
        time.append(dcc[k].sleep)
    dcc = np.concatenate(dcc); rat = relable(np.concatenate(rat)); session = np.concatenate(session)
    isa = relable(np.concatenate(isa)) - 0.5; time = np.concatenate(time) - 0.5

    df = pd.DataFrame({"y": dcc, "rat": rat, "session": session, "isa": isa, "time": time})
    df["sess_id"] = df["subject"].astype(str) + "_" + df["session"].astype(str)
    m = smf.mixedlm("y ~ cond * time", df,
                    groups="subject",                          # random intercept per subject
                    vc_formula={"session": "0 + C(sess_id)"}   # session nested in subject
                    ).fit(reml=True)
    print(m.summary())


        # isau.pair_violin(tuple(dcc), tuple(d.rat for d in dcc), ax=ax[2*i+j], color=isau.paperColors((r,'shuffle')),
        #         label=('ISA      ','      nISA'), ylim=(0,0.5), ylim2=(0,0.6))
        # ax[2*i+j].set_title('post-task' if j else 'baseline', size=6)
        # if j != 0:
        #     ax[2*i+j].spines.left.set_visible(False); ax[2*i+j].set(yticks=())
        # if j != 1: ax[2*i+j].set(xticks=())
ax[0].set(ylabel='DCC')



fig.show()
fma.plotting.saveFigure(fig, output_file, ('png','svg','pdf'))
isau.save_axes_figs(axs, fig, output_file)