import fmatoolbox as fma
import isautilities as isau
import numpy as np
import xarray as xr
import pandas as pd
import statsmodels.formula.api as smf
import pathlib
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig5'
results = isau.load_results(res_root / 'criticalExponents.pkl')
FW = 17.5; FH = 10.9; LEFT = 1; TOP = 0.5; GAP = 1
fig = fma.plotting.figure(n=0, size=(FW,FH))

# PFC criticality
fit1 = results['fit_ISA']; exp1 = results['exp_ISA']; fit2 = results['fit_ISA_s2']; exp2 = results['exp_ISA_s2']
base = fit2.base; states = fit2.state.values; regs = fit2.reg.values
fit = xr.concat((fit1,fit2), dim='sleep', join='outer'); fit = fit.sel(base=base)
exp = xr.concat((exp1,exp2), dim='sleep', join='outer'); exp = exp.sel(base=base)
relable = lambda x : np.unique(x, return_inverse=True)[1]
ax = []
for i, r in enumerate(regs):
    dcc = []; rat = []; session = []; isa = []; time = []; valid = []
    for sl in range(2):
        for k, s in enumerate(states):
            h = fit.sel(sleep=sl, reg=r, state=s, val='p').values < 0.05 # failed sessions get h False since p is nan
            signif = np.all(h, axis=1) # sessions where the fit was performed and significant for S and D
            dcc.append(exp.sel(sleep=sl, reg=r, state=s, val='DCC').copy(deep=True))
            dcc[-1][~signif] = np.nan
            # save metadata
            rat.append(dcc[k].rat.values)
            session.append(dcc[k].base)
            isa.append(np.full(dcc[k].shape, k))
            time.append(np.full(dcc[k].shape, sl))
            valid.append(signif.astype(int))
    valid = np.mean(valid, axis=0) > 0.5 # sessions where at least half conditions were significant

    for k in range(len(dcc)):
        dcc[k] = dcc[k][valid]; rat[k] = rat[k][valid]; session[k] = session[k][valid]; isa[k] = isa[k][valid]; time[k] = time[k][valid]
        dcc[k][np.isnan(dcc[k])] = dcc[k][~np.isnan(dcc[k])].max() # input non significant fits with highest observed distance

    ax.append(isau.axes(fig, LEFT+i*5.3, TOP, 5,2.5))
    fma.plotting.boxPlot(dcc, mode='scatter', ax=ax[-1], color=isau.paperColors(r))
    #isau.pair_violin((dcc[0],dcc[1]), (rat[0],rat[1]), ax=ax[-1], color=isau.paperColors((r,'shuffle')), label=('ISA','nISA'))
    #isau.pair_violin((dcc[2],dcc[3]), (rat[2],rat[3]), center=2, ax=ax[-1], color=isau.paperColors((r,'shuffle')), label=('ISA','nISA'))
    ax[-1].set(xlim=(-1,4))

    # relabel rats and sessions in {0, N}, map condition and time to +/- 0.5
    dcc = np.concatenate(dcc); rat = relable(np.concatenate(rat)); session = relable(np.concatenate(session))
    isa = np.concatenate(isa) - .5; time = np.concatenate(time) - .5
    df = pd.DataFrame({"dcc": dcc, "rat": rat, "session": session, "isa": isa, "time": time})
    m = smf.mixedlm("dcc ~ isa * time", df,
                    groups="rat", re_formula="1",              # random intercept per rat
                    vc_formula={"session": "0 + C(session)"}   # session nested in subject
                    ).fit(reml=True)
    print(r.upper(), m.summary())


fig.show()