import fmatoolbox as fma
import isautilities as isau
import pathlib
import numpy as np
import scipy.stats

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'Fig4'
LEFT = 0.8; TOP = 0.5; VSPACE = 1.2
fig = fma.plotting.figure(n=0, size=(17.5,11.5))


# criticality in ISA
results = isau.load_results(res_root / 'criticalExponents.pkl')
width = 1.5
ax = [isau.axes(fig, LEFT+i*(width+0.5), TOP, width, 2.5) for i in range(2)]
fit = results['fit_ISA']; exponents = results['exp_ISA']
regs = ('pfc','hpc'); states = fit.state.values
axs = {f'dcc_{regs[i]}': a for i, a in enumerate(ax)}
dcc = [[]] * 2
ylim = (0,0.5); yticks = (0,0.2,0.4,0.6)
for i, r in enumerate(regs):
    for j, s in enumerate(states):
        valid = np.all(~np.isnan(fit.sel(reg=r,state=s,val='p')), axis=1) # sessions where a fit was performed
        h = fit.sel(reg=r,state=s,val='p')[valid].values < 0.05
        signif = np.all(h, axis=1) # sessions where the fit was successful
        dcc[j] = exponents.sel(reg=r,state=s,val='DCC')[valid]
        dcc[j][~signif] = np.nan
    valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) # sessions where at least one condition was significant
    for j, s in enumerate(states):
        dcc[j] = dcc[j][valid]
        dcc[j][np.isnan(dcc[j])] = dcc[j][~np.isnan(dcc[j])].max() # input non significant fits with highest observed distance

    isau.pair_violin((dcc[0],dcc[1]), (dcc[0].rat,dcc[1].rat), color=isau.paperColors((r,'shuffle')), label=('ISA','nISA'),
                     ax=ax[i], ylim=ylim, yticks=yticks)
    ax[i].set(ylim=(0,0.6))
ax[0].set(ylabel='DCC')


# xCorr    SHOW THAT PEAKS ARE ~ZERO FOR nISA AND >0 FOR ISA ??
results = isau.load_results(res_root / 'unitLagCorr.pkl')
left, top, width, height = isau.pos(ax[-1], fig)
ax = isau.axes(fig, left+width+VSPACE, top, 2.8, height)
axs['corr'] = ax
corr = results['corr']
reg = 'hpc-pfc'; color = isau.paperColors(('shuffle','nr'))
ax.axvline(0, ls='--', color='k', lw=0.5, zorder=0)
for j, c in enumerate(['nisa','isa']):
    reg_data = corr.sel(cond=c, reg=reg)
    reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True)
    fma.plotting.semPlot(reg_data.lag, reg_data, color=color[j])
ax.set(xlabel='PFC to HPC lag (s)')
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(5)
# stats
for j, c in enumerate(['nisa','isa']):
    reg_data = corr.sel(cond=c, reg=reg)
    reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True)
    p = scipy.stats.wilcoxon(reg_data, axis=0, alternative='greater', nan_policy='omit').pvalue
    _, p = fma.analysis.holmBonferroni(p, return_reject=True)
    fma.plotting.pHorzLine(p, reg_data.lag, color=color[j], lw=2)
ax.set(xticks=(-0.15,-0.05,0.05,0.15), ylabel=r'$\langle$Δ Spearman corr$\rangle$', yticks=np.array((-0.3,0,0.3,0.6,0.9))*1e-3)


def plotCMap(x,ax=None,reg=None,**kwargs):
    if ax is not None and reg is not None:
        reg = str(reg)
        ax.set(xlabel=f'{reg.split('-')[1].upper()} lag (s)',ylabel=f'{reg.split('-')[0].upper()} lag (s)')
        out = fma.plotting.plotColorMap(x, x=x.lag, y=x.target, ax=ax, interpolation='none', **kwargs)
        if isinstance(out, tuple):
            out[1].ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
            out[1].ax.yaxis.get_offset_text().set_fontsize(5)
    return out


# ripples xCorr
results = isau.load_results(res_root / 'unitPECorr.pkl')
left, top, width, height = isau.pos(ax, fig)
ax = [isau.axes(fig, left+width+VSPACE+i*(height+VSPACE), top, height, height) for i in range(2)]
for n, a in zip(('xcorr_rip_ISA','xcorr_rip_nISA'), ax): axs[n] = a
corr = results['rip']
vmin = 0.0010; vmax = 0.0033
label = [None,r'$\langle$Δ Spearman corr$\rangle$']
reg = corr.reg.values[1]
for a, isa in enumerate([True,False]):
    reg_data = corr.sel(reg=reg, isa=isa)
    reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True)
    # significance
    w = scipy.stats.wilcoxon(reg_data, axis=0, nan_policy='omit', alternative='greater').pvalue
    _, h = fma.analysis.holmBonferroni(w, return_reject=True)
    # plot result
    out = plotCMap(reg_data.mean(dim='rat',skipna=True).T, vmin=vmin, vmax=vmax, smooth=.2, alpha=np.where(h,1.0,0.5),
                   bar=label[a], ax=ax[a], reg=reg, aspect='preserve')
    ax[a].set_title('ISA' if isa else 'nISA', fontsize=7, pad=2)
out[1].ax.set_box_aspect(None)
bb = ax[-1].get_position()
dy = 0.02
out[1].ax.set_position([bb.xmax+0.03, bb.ymin+dy, 0.01, bb.height-2*dy])
axs['xcorr_rip_cb'] = out[1].ax

# bootstrap upper-right quadrant
left, top, width, height = isau.pos(ax[-2], fig)
ax = isau.axes(fig, left, top+height+VSPACE, 2.8, height)
def mask_quadrant(x,a=-2,b=-2,n=7):
    i, j = np.indices(x.shape[-2:])
    mask = (i >= len(i)//2 + a) & (j >= len(j)//2 + b) & (i - j > -n) & (i - j < n)
    x = np.where(mask,x,np.nan) # upper-right quadrant
    return x
valid = corr.rat != 11
reg_data = corr.sel(reg=corr.reg.values[1],rat=valid)
reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True) # increase from shuffle
reg_data = np.nanmax(mask_quadrant(reg_data.values),(2,3)) # upper-right quadrant
isau.pair_violin((reg_data[:,0],reg_data[:,1]), (corr.rat.sel(rat=valid),)*2, ax, mode='diff', color=isau.paperColors(('ripples','shuffle')),
                     label=('ISA','nISA'), ylim=(0,0.01), yticks=(0,0.005,0.01), ylabel=r'$\langle$Δ Spearman corr$\rangle$')
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
ax.yaxis.get_offset_text().set_fontsize(5)
axs['peak_rip'] = ax


# deltas xCorr
left, top, width, height = isau.pos(ax, fig)
ax = [isau.axes(fig, LEFT+i*(height+VSPACE), top, height, height) for i in range(2)]
for n, a in zip(('xcorr_del_ISA','xcorr_del_nISA'), ax): axs[n] = a
corr = results['delta']
vmin = 0.0005; vmax = 0.0027
label = [None,r'$\langle$Δ Spearman corr$\rangle$']
reg = corr.reg.values[1]
for a, isa in enumerate([True,False]):
    reg_data = corr.sel(reg=reg, isa=isa)
    reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True)
    # significance
    w = scipy.stats.wilcoxon(reg_data, axis=0, nan_policy='omit', alternative='greater').pvalue
    _, h = fma.analysis.holmBonferroni(w, return_reject=True)
    # plot result
    out = plotCMap(reg_data.mean(dim='rat',skipna=True).T, vmin=vmin, vmax=vmax, smooth=.2, alpha=np.where(h,1.0,0.5),
                   bar=label[a], ax=ax[a], reg=reg, aspect='preserve', cmap='plasma')
    ax[a].set_title('ISA' if isa else 'nISA', fontsize=7, pad=2)
out[1].ax.set_box_aspect(None)
bb = ax[-1].get_position()
dy = 0.02
out[1].ax.set_position([bb.xmax+0.03, bb.ymin+dy, 0.01, bb.height-2*dy])
axs['xcorr_del_cb'] = out[1].ax






# spindles xCorr
_, top, width, height = isau.pos(ax[-1], fig)
ax = [isau.axes(fig, LEFT+i*(height+VSPACE), top+height+VSPACE, height, height) for i in range(2)]
for n, a in zip(('xcorr_spin_ISA','xcorr_spin_nISA'), ax): axs[n] = a
corr = results['spin']
vmin = 0.000; vmax = 0.002
label = [None,r'$\langle$Δ Spearman corr$\rangle$']
reg = corr.reg.values[1]
for a, isa in enumerate([True,False]):
    reg_data = corr.sel(reg=reg, isa=isa)
    reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True)
    # significance
    w = scipy.stats.wilcoxon(reg_data, axis=0, nan_policy='omit', alternative='greater').pvalue
    _, h = fma.analysis.holmBonferroni(w, return_reject=True)
    # plot result
    out = plotCMap(reg_data.mean(dim='rat',skipna=True).T, vmin=vmin, vmax=vmax, smooth=.2, alpha=np.where(h,1.0,0.5),
                   bar=label[a], ax=ax[a], reg=reg, aspect='preserve', cmap='cividis')
    ax[a].set_title('ISA' if isa else 'nISA', fontsize=7, pad=2)
out[1].ax.set_box_aspect(None)
bb = ax[-1].get_position()
dy = 0.02
out[1].ax.set_position([bb.xmax+0.03, bb.ymin+dy, 0.01, bb.height-2*dy])
axs['xcorr_spin_cb'] = out[1].ax




fig.show()
fma.plotting.saveFigure(fig, output_file, ['png','svg'])
isau.save_axes_figs(axs, fig, output_file)