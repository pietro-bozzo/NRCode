import fmatoolbox as fma
import isautilities as isau
import pathlib
import numpy as np
import scipy.stats

def mask_quadrant(x, a=-2, b=-2, c=15, d=15, n=5):
    i, j = np.indices(x.shape[-2:])
    mask = (i <= len(i)//2 + c) & (j <= len(j)//2 + d) & (i >= len(i)//2 + a) & (j >= len(j)//2 + b) & (i - j > -n) & (i - j < n)
    x = np.where(mask, x, np.nan) # upper-right quadrant
    return x

def main():

    # parameters
    froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
    res_root = froot / 'Data'
    output_file = froot / 'FinalFigures' / 'Fig4'
    LEFT = 0.8; TOP = 0.5; VSPACE = 1.2; SMOOTH = 0.3; SAVE = True
    fig = fma.plotting.figure(n=0, size=(17.5,7.5))


    # criticality in ISA
    results = isau.load_results(res_root / 'criticalExponents.pkl')
    width = 1
    ax = [isau.axes(fig, LEFT+i*width, TOP, width, 2.5) for i in range(2)]
    fit = results['fit_ISA']; exponents = results['exp_ISA']
    regs = ('pfc','hpc'); states = fit.state.values
    axs = {f'dcc_{regs[i]}': a for i, a in enumerate(ax)}
    dcc = [[]] * 2
    ylim = (0,0.5); yticks = (0,0.2,0.4,0.6)
    for i, r in enumerate(regs):
        for j, s in enumerate(states):
            h = fit.sel(reg=r, state=s, val='p').values < 0.05 # failed sessions get h False since p is nan
            signif = np.all(h, axis=1) # sessions where the fit was performed and significant for S and D
            dcc[j] = exponents.sel(reg=r, state=s, val='DCC').copy(deep=True)
            dcc[j][~signif] = np.nan
        valid = ~np.isnan(dcc[0]) | ~np.isnan(dcc[1]) # keep sessions where at least one condition was significant
        for j, s in enumerate(states):
            dcc[j] = dcc[j][valid]
            dcc[j][np.isnan(dcc[j])] = dcc[j][~np.isnan(dcc[j])].max() # input non significant fits with highest observed distance

        isau.pair_violin((dcc[0],dcc[1]), (dcc[0].rat,dcc[1].rat), color=isau.paperColors((r,'shuffle')), label=('ISA','nISA'),
                         alternative='greater',ax=ax[i], ylim=ylim, yticks=yticks)
        ax[i].set(ylim=(0,0.6))
    ax[0].set(ylabel='DCC'); ax[1].set(yticks=());  ax[1].spines.left.set_visible(False)


    # xCorr
    results = isau.load_results(res_root / 'unitLagCorr.pkl')
    left, top, width, height = isau.pos(ax[-1])
    ax = isau.axes(fig, left+width+VSPACE, top, 2.8, height)
    axs['corr'] = ax
    corr = results['corr_s1']; reg = 'hpc-pfc'
    color = isau.paperColors(('shuffle','nr'))
    ax.axvline(0, ls='--', color='k', lw=0.5, zorder=0)
    for j, cond in enumerate(['nisa','isa']):
        c = corr.sel(cond=cond, reg=reg)
        c = c.sel(shuffle=False) - c.sel(shuffle=True).mean(dim='shuffle')
        fma.plotting.semPlot(c.lag, c, smooth=0.7, color=color[j])
    ax.set(xlabel='PFC to HPC lag (s)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax.yaxis.get_offset_text().set_fontsize(5)
    # stats
    p = []
    for j, cond in enumerate(['nisa','isa']):
        c = corr.sel(cond=cond, reg=reg)
        p.append(fma.statistics.maxStatisticTest(c.sel(shuffle=False), c.sel(shuffle=True).transpose(..., 'shuffle'),
                                                 group=c.rat))
    fma.plotting.pHorzLine(np.array(p).T, c.lag, color=color, lw=1)
    ax.set(xticks=(-0.15,-0.05,0.05,0.15), ylabel=r'$\langle$Δ Spearman corr$\rangle$',
           yticks=np.array((-0.2,0,0.2,0.4,0.6,0.8))*1e-3)

    # peak corr time
    ax = isau.axes(fig, 0.4, 0, 0.85, 0.5, ax, 'lh')
    axs['corr_peak'] = ax
    peaks = []
    for j, cond in enumerate(['isa','nisa']):
        c = corr.sel(cond=cond, reg=reg)
        c = c.sel(shuffle=False) - c.sel(shuffle=True).mean(dim='shuffle')
        peak = c.dropna(dim='rat', how='all').argmax(dim='lag')
        peaks.append(corr.lag[peak])
    isau.pair_violin(tuple(peaks), (peaks[0].rat,peaks[1].rat), color=isau.paperColors(('nr','shuffle')), label=('ISA','nISA'),
                     alternative='less', ylim2=(-0.05,0.05), ax=ax)
    #isau.bayesHier(tuple(peaks), (peaks[0].rat,peaks[1].rat))


    # ripples xCorr
    results = isau.load_results(res_root / 'unitPECorr.pkl')
    left, top, width, height = isau.pos(axs['corr'])
    ax = [isau.axes(fig, left+width+VSPACE*2+i*(height+0.7), top, height, height) for i in range(2)]
    for n, a in zip(('xcorr_rip_ISA','xcorr_rip_nISA'), ax): axs[n] = a
    corr = results['rip']; corr = corr.sel(shuffle=False) - corr.sel(shuffle=True).mean(dim='shuffle') # increase from shuffle
    vmin = 0.001; vmax = 0.0027
    label = [None,r'$\langle$Δ Spearman corr$\rangle$']
    reg = 'hpc-pfc'
    AA = []
    for a, isa in enumerate([True,False]):
        reg_data = corr.sel(reg=reg, isa=isa)
        #AA.append( mask_quadrant(reg_data.mean(dim='rat',skipna=True)) )
        AA.append(reg_data.mean(dim='rat',skipna=True))
        out = isau.plotCMap(AA[-1].T, vmin=vmin, vmax=vmax, smooth=SMOOTH, bar=label[a], ax=ax[a], reg=reg, aspect='preserve')
        ax[a].set_title('ISA' if isa else 'nISA', fontsize=7, pad=2)
    ax[1].set(ylabel='')
    out[1].ax.set_box_aspect(None)
    bb = ax[-1].get_position()
    dy = 0.1
    out[1].ax.set_position([bb.xmax+0.01, bb.ymin, 0.01, bb.height-dy])
    axs['xcorr_rip_cb'] = out[1].ax

    # bootstrap upper-right quadrant
    ax = isau.axes(fig, 1.5, 0, 0.85, 0.5, ax[-1], 'lh')
    valid = corr.rat != 11
    reg_data = corr.sel(reg=reg, rat=valid)
    reg_data = np.nanmax(mask_quadrant(reg_data),(2,3)) # upper-right quadrant
    #reg_data = np.nanmean(reg_data,(2,3))
    isau.pair_violin((reg_data[:,0],reg_data[:,1]), (corr.rat.sel(rat=valid),)*2, ax, mode='diff', color=isau.paperColors(('ripples','shuffle')),
                         label=('ISA','nISA'), ylim=(0,0.01), yticks=(0,0.005,0.01), ylabel=r'$\langle$Δ Spearman corr$\rangle$')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax.yaxis.get_offset_text().set_fontsize(5)
    axs['peak_rip'] = ax



    # deltas xCorr
    left, top, width, height = isau.pos(axs['xcorr_rip_ISA'])
    ax = [isau.axes(fig, LEFT+i*(height+VSPACE), top+height+VSPACE, height, height) for i in range(2)]
    for n, a in zip(('xcorr_del_ISA','xcorr_del_nISA'), ax): axs[n] = a
    corr = results['delta']
    vmin = 0.0005; vmax = 0.0027
    label = [None,r'$\langle$Δ Spearman corr$\rangle$']
    for a, isa in enumerate([True,False]):
        reg_data = corr.sel(reg=reg, isa=isa)
        reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True).mean(dim='shuffle')
        out = isau.plotCMap(reg_data.mean(dim='rat',skipna=True).T, vmin=vmin, vmax=vmax, smooth=SMOOTH,
                       bar=label[a], ax=ax[a], reg=reg, aspect='preserve', cmap='plasma')
        ax[a].set_title('ISA' if isa else 'nISA', fontsize=7, pad=2)
    out[1].ax.set_box_aspect(None)
    bb = ax[-1].get_position()
    dy = 0.02
    out[1].ax.set_position([bb.xmax+0.03, bb.ymin+dy, 0.01, bb.height-2*dy])
    axs['xcorr_del_cb'] = out[1].ax


    # spindles xCorr
    left, top, width, height = isau.pos(axs['xcorr_rip_ISA'])
    ax = [isau.axes(fig, left+i*(height+VSPACE), top+height+VSPACE, height, height) for i in range(2)]
    for n, a in zip(('xcorr_spin_ISA','xcorr_spin_nISA'), ax): axs[n] = a
    corr = results['spin']
    vmin = 0.0001; vmax = 0.00055
    label = [None,r'$\langle$Δ Spearman corr$\rangle$']
    reg = corr.reg.values[1]
    for a, isa in enumerate([True,False]):
        reg_data = corr.sel(reg=reg, isa=isa)
        reg_data = reg_data.sel(shuffle=False) - reg_data.sel(shuffle=True).mean(dim='shuffle')
        out = isau.plotCMap(reg_data.mean(dim='rat',skipna=True).T, vmin=vmin, vmax=vmax, smooth=SMOOTH,
                       bar=label[a], ax=ax[a], reg=reg, aspect='preserve', cmap='cividis')
        ax[a].set_title('ISA' if isa else 'nISA', fontsize=7, pad=2)
    out[1].ax.set_box_aspect(None)
    bb = ax[-1].get_position()
    dy = 0.02
    out[1].ax.set_position([bb.xmax+0.03, bb.ymin+dy, 0.01, bb.height-2*dy])
    axs['xcorr_spin_cb'] = out[1].ax



    fig.show()
    SAVE and fma.plotting.saveFigure(fig, output_file, ('png','svg','pdf'))
    SAVE and isau.save_axes_figs(axs, fig, output_file)

if __name__ == "__main__":
    main()