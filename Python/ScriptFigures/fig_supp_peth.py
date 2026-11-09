import fmatoolbox as fma
import isautilities as isau
import pathlib

# parameters
froot = pathlib.Path().cwd().parent.parent.parent / 'Results'
res_root = froot / 'Data'
output_file = froot / 'FinalFigures' / 'FigSuppPETH'
LEFT = 1; TOP = 0.5; HSPACE = 0.7; VSPACE = 1
axs = {}
fig = fma.plotting.figure(n=0, size=(8.5,9.5))

# PETH OFF-ON
results = isau.load_results(res_root / 'HpcPfcPETH.pkl')
width = 2.9; height = 2
ax = [isau.axes(fig, LEFT, TOP+i*(VSPACE+height), width, height) for i in range(3)]
for n, a in zip(('peth_rip_off_on','peth_delta_off_on','peth_spin_off_on'), ax): axs[n] = a
peth = results['peth']
ylim = ((0.008,0.03), (-0.0002,0.015), (-0.0002,0.02))
yticks = ((0.01,0.02,0.03), (0,0.005,0.010,0.015), (0,0.01,0.02))
for i, e in enumerate(peth.ev.values):
    ax[i].axvline(0, ls='--', lw=0.5, color='gray')
    tc = peth.sel(ev=e, tr='off-on')
    fma.plotting.semPlot(tc.t, tc, color=isau.paperColors(e), ax=ax[i], smooth=1.5)
    ax[i].set(xlabel='time from OFF-ON (s)', ylabel='rate (Hz)', ylim=ylim[i], yticks=yticks[i])
    ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax[i].yaxis.get_offset_text().set_fontsize(5)



ax = [isau.axes(fig, HSPACE, i*(VSPACE+height), 1, 1, ax[0], 'lwh') for i in range(3)]
for n, a in zip(('peth_rip_on_off','peth_delta_on_off','peth_spin_on_off'), ax): axs[n] = a
for i, e in enumerate(peth.ev.values):
    ax[i].axvline(0, ls='--', lw=0.5, color='gray')
    tc = peth.sel(ev=e, tr='on-off')
    fma.plotting.semPlot(tc.t, tc, color=isau.paperColors(e), ax=ax[i], smooth=1.5)
    ax[i].set(xlabel='time from OFF-ON (s)', ylim=ylim[i], yticks=yticks[i])
    ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    ax[i].yaxis.get_offset_text().set_fontsize(5)



fig.show()
fma.plotting.saveFigure(fig, output_file, ('png','svg','pdf'))
isau.save_axes_figs(axs, fig, output_file)