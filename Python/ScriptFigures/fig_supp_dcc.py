import numpy as np
import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'FigSuppDCC'
LEFT = 0.7; TOP = 0.7
axs = {}
fig = fma.plotting.figure(n=0, size=(12,5))
results = isau.load_results(res_root / 'criticalExponents.pkl')


# DCC across brain states
fit = results['fit_states']; exponents = results['exp_states']; states = fit.state.values; regs = fit.reg.values
ax = [isau.axes(fig, LEFT+i*(3+1), TOP,3,2) for i in range(len(regs))]
axs['dcc'] = ax
for i, r in enumerate(regs):
    dcc = []
    for j, s in enumerate(states):
        valid = np.all(~np.isnan(fit.sel(reg=r,state=s,val='p')), axis=1) # sessions where a fit was performed
        h = fit.sel(reg=r,state=s,val='p')[valid].values < 0.05
        signif = np.all(h, axis=1) # sessions where the fit was successful
        dcc.append(exponents.sel(reg=r,state=s,val='DCC'))
        dcc[-1][~valid] = np.nan
        dcc[-1][valid][~signif] = np.nan
    valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) | ~np.isnan(dcc[2]) # sessions where at least one condition was significant
    for j in range(len(dcc)):
        dcc[j] = dcc[j][valid]
        dcc[j][np.isnan(dcc[j])] = dcc[j][~np.isnan(dcc[j])].max() # input non significant fits with highest observed distance
    b = fma.statistics.hierarchicalBootstrap(tuple(dcc), tuple(d.rat for d in dcc))
    fma.plotting.boxPlot(b, mode='violin', color=isau.paperColors(states), label=('nREM','REM','wake'), ax=ax[i])
    p = fma.statistics.hierarchicalBootstrapTest(tuple(dcc), tuple(d.rat for d in dcc), paired=True)
    ax[i].set(ylim=(0,0.5))
    fma.plotting.pBar(p[['Condition 1', 'Condition 2', 'Corrected p-value']], ax=ax[i])
    ax[i].set(xlim=(-0.5,2.5), ylim=(0,0.6))
    ax[i].set_title(r.upper(), size=6)
ax[0].set(ylabel='DCC')


fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])