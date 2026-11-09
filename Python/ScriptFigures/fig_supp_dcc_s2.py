import numpy as np
import fmatoolbox as fma
import isautilities as isau
import xarray as xr
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'FigSuppStatesS2'
LEFT = 0.9; TOP = 0.5; GAP = 1.
axs = {}
fig = fma.plotting.figure(n=0, size=(17.5,8))
results = isau.load_results(res_root / 'criticalExponents.pkl')


# DCC in ISA vs nISA     PFC TO SUPP!!
fit1 = results['fit_ISA']; exponents1 = results['exp_ISA']; fit2 = results['fit_ISA_s2']; exponents2 = results['exp_ISA_s2']
print(fit2.base)
states = fit1.state.values; regs = fit1.reg.values
width = 1.5
ax = []; left = LEFT
for i in range(len(regs) * len(states)):
    ax.append(isau.axes(fig, left, TOP, width,2))
    left += width + GAP * int(i % 2)
axs['dcc_isa_s2'] = ax
for i, r in enumerate(regs):
    for j, (fit, exp) in enumerate(zip((fit1,fit2), (exponents1,exponents2))):
        dcc = []
        for k, s in enumerate(states):
            h = fit.sel(reg=r, state=s, val='p').values < 0.05 # failed sessions get h False since p is nan
            signif = np.all(h, axis=1) # sessions where the fit was performed and significant for S and D
            dcc.append(exp.sel(reg=r, state=s, val='DCC').copy(deep=True))
            dcc[-1][~signif] = np.nan
        valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) # sessions where at least one condition was significant
        for k in range(len(dcc)):
            dcc[k] = dcc[k][valid]
            dcc[k][np.isnan(dcc[k])] = dcc[k][~np.isnan(dcc[k])].max() # input non significant fits with highest observed distance

        isau.pair_violin(tuple(dcc), tuple(d.rat for d in dcc), ax=ax[2*i+j], color=isau.paperColors((r,'shuffle')),
                label=('ISA      ','      nISA'), ylim=(0,0.5), ylim2=(0,0.6))
        ax[2*i+j].set_title('post-task' if j else 'baseline', size=6)
        if j != 0:
            ax[2*i+j].spines.left.set_visible(False); ax[2*i+j].set(yticks=())
        if j != 1: ax[2*i+j].set(xticks=())
ax[0].set(ylabel='DCC')


# DCC across brain states
fit1 = results['fit_states']; exponents1 = results['exp_states']; fit2 = results['fit_states_s2']; exponents2 = results['exp_states_s2']
fit = xr.concat((fit1,fit2), dim='sleep', join='inner'); exp = xr.concat((exponents1,exponents2), dim='sleep', join='inner')
states = fit1.state.values; regs = fit1.reg.values; state_names = ('nREM','REM','wake')
width = 1.5
ax = []; left = LEFT
for i in range(len(regs) * len(states)):
    ax.append(isau.axes(fig, left, TOP+3.5, width,2))
    left += width + GAP * int(i%3==2)
axs['dcc_states_s2'] = ax
for i, r in enumerate(regs):
    for k, s in enumerate(states):
        dcc = []
        for sl in range(2):
            h = fit.sel(sleep=sl, reg=r, state=s, val='p').values < 0.05 # failed sessions get h False since p is nan
            signif = np.all(h, axis=1) # sessions where the fit was performed and significant for S and D
            dcc.append(exp.sel(sleep=sl, reg=r, state=s, val='DCC').copy(deep=True))
            dcc[-1][~signif] = np.nan
        valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) # sessions where at least one condition was significant
        for j in range(len(dcc)):
            dcc[j] = dcc[j][valid]
            dcc[j][np.isnan(dcc[j])] = dcc[j][~np.isnan(dcc[j])].max() # input non significant fits with highest observed distance

        isau.pair_violin(tuple(dcc), tuple(d.rat for d in dcc), ax=ax[3*i+k], color=isau.paperColors((r,'shuffle')),
                label=('baseline      ','      post-task'), ylim=(0,0.5), ylim2=(0,0.6))
        ax[3*i+k].set_title(f'{r.upper()} {state_names[k]}', size=6)
        if k != 0:
            ax[3*i+k].spines.left.set_visible(False); ax[3*i+k].set(yticks=())
        if k != 1: ax[3*i+k].set(xticks=())
ax[0].set(ylabel='DCC')


fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])