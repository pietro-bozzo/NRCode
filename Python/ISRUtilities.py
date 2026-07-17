import numpy as np
import fmatoolbox as fma
import matplotlib.axes as mpla
from matplotlib.pyplot import yticks


def paperColors(i,alpha='ff'):
    # get colors for Nucleus Reuniens Infra-Slow Rhythm paper
    #
    # arguments:
    #     i        (n) list of int, color indeces, interpreted modulo[number of colors]
    #     alpha    string = 'ff', transparency
    #
    # output:
    #     colors    (n) list of string, hex colors

    # IBM colorblind safe palette
    #colors = ["#648fff","#ffb000","#dc267f","#785ef0","#fe6100"] * 3

    colors = {'hpc':           '#87b69e', # 1. HPC
              'nr':            '#a076b9', # 2. NR
              'pfc':           "#fc7860", # 3. PFC
              'th':            "#f4ce32", # 4. TH
              'sws':           "#ecc8e3", # 5. nREM
              'rem':           "#a8c7ec", # 6. REM
              'wake':          "#fff3b6", # 7. wake
              'on':            "#d8e9bc", # 8. NR "up state"
              'off':           "#bee5eb", # 9. NR "down state"
              'ripples':       "#1f99a3", # 10. ripples
              'deltas':        "#c80000", # 11. deltas
              'spindles':      "#2a619e", # 12. spindles
              'shuffle':       "#888888", # 13. shuffle
              'v1':            "#f4ce32", # 4bis. V1
              'amy':           "#f4ce32", # 4bis. AMY
              'lamy':          "#f4ce32", # 4bis. left AMY
              'ramy':          "#f4ce32", # 4bis. right AMY
              'lvmpfc':        "#f4ce32", # 4bis. left ventral PFC
              'rvmpfc':        "#f4ce32", # 4bis. right ventral PFC
              'ldmpfc':        "#f4ce32", # 4bis. left dorsal PFC
              'rdmpfc':        "#f4ce32", # 4bis. right dorsal PFC
              'lvhpc':         "#f4ce32", # 4bis. left ventral HPC
              'rvhpc':         "#f4ce32", # 4bis. right ventral HPC
              'ldhpc':         "#f4ce32", # 4bis. left dorsal HPC
              'rdhpc':         "#f4ce32", # 4bis. right dorsal HPC
              'deltawaves':    "#c80000", # 11bis. deltas
              }
    n_colors = len(colors)
    for j, key in enumerate(tuple(colors)):
        colors[j] = colors[key]

    if isinstance(i,int):
        return colors[i % n_colors] + alpha
    elif isinstance(i,str):
        out = colors[i.lower()] if i.lower() in colors else colors['th']
        return out + alpha
    else:
        out = []
        for j in i:
            if isinstance(j,int):
                out.append(colors[j % n_colors] + alpha)
            else:
                o = colors[j.lower()] if j.lower() in colors else colors['th']
                out.append(o + alpha)
        return out


def ISAxes(axs, double:bool=None, **ax_kwargs):

    if isinstance(axs,mpla._axes.Axes):
        axs = [axs]
    elif isinstance(axs,np.ndarray):
        axs = axs.ravel()
    if double is None: double = False

    on_color = '#444444'
    off_color = '#c9a6f7'
    lw = 3

    for ax in axs:
        if ax.name == 'polar':
            ax.set_theta_zero_location('S')
            ax.spines['polar'].set_visible(False)
            ax.set_xticks([0,np.pi/2,4*np.pi/6,np.pi,3*np.pi/2,10*np.pi/6],['','','ON','','','OFF'],fontsize=12)
            ax.set(yticks=[])

            labels = ax.get_xticklabels()
            labels[2].set_color(on_color)
            labels[5].set_color(off_color)

            gridlines = ax.xaxis.get_gridlines()
            for i in [0,1,3,4]:
                gridlines[i].set_zorder(i)
            for i in [2,5]:
                gridlines[i].set_visible(False)

            ax.set(**ax_kwargs)
            ylim = ax.get_ylim()

            theta = np.linspace(0,np.pi,200)
            ax.plot(theta,np.full_like(theta,ylim[1]),lw=lw,color=on_color,clip_on=False,zorder=5)
            ax.plot(theta+np.pi,np.full_like(theta,ylim[1]),lw=lw,color=off_color,clip_on=False,zorder=6)
            ax.scatter(np.pi,ylim[1],s=80,color=on_color,marker='<',clip_on=False,zorder=7)
            ax.scatter(0,ylim[1],s=80,color=off_color,marker='>',clip_on=False,zorder=8)
            ax.set_ylim(ylim)

        else:
            ax.spines['bottom'].set_visible(False)

            x_ticks = np.arange(0,(double+1)*2*np.pi,np.pi) + np.pi/2
            ax.set_xticks(x_ticks,['ON','OFF']*(double+1),fontsize=12)
            labels = ax.get_xticklabels()
            for i in range(double+1):
                labels[2*i].set_color(on_color)
                labels[2*i+1].set_color(off_color)
            ax.tick_params(axis='x',length=0)

            ax.set(**ax_kwargs)
            ylim = ax.get_ylim()
            x_coord = [0,np.pi]
            if double:
                x_coord.extend([np.nan,2*np.pi,3*np.pi])
            ax.plot(x_coord,np.full_like(x_coord,ylim[0]),lw=lw,color=on_color,clip_on=False)
            ax.plot(np.array(x_coord)+np.pi,np.full_like(x_coord,ylim[0]),lw=lw,color=off_color,clip_on=False)
            ax.set_ylim(ylim)

    return


def ISCycles(isa,on,return_off=False):

    isa = np.array(isa,ndmin=2)
    on = np.array(on,ndmin=2)
    if isa.size == 0 or on.size == 0:
        return np.zeros((0,2)), np.zeros((0,2)) if return_off else np.zeros((0,2))

    _, is_inside = fma.general.restrict(on[1:,0] - .0001, isa,s_ind=True) # is_inside[i] is False if on[i+1,0] cannot be used as end of a cycle (or OFF)
    is_inside = np.append(is_inside,False) # last element of cycles[:,1] (or off[:,1]) is always to change

    cycles = np.stack((on[:,0], np.roll(on[:,0],-1)), axis=1)
    isa_idx = np.searchsorted(isa[:,1], cycles[:,1][~is_inside]) - 1  # find which ISA end will replace invalid ON ends
    cycles[:,1][~is_inside] = isa[isa_idx,1]

    if return_off:
        off = np.stack((on[:,1], np.roll(on[:,0],-1)), axis=1)
        isa_idx = np.searchsorted(isa[:,1], off[:,1][~is_inside]) - 1  # find which ISA end will replace invalid ON ends
        off[:,1][~is_inside] = isa[isa_idx,1]
        return cycles, off

    return cycles


def ISTransitions(isa,on,min_duration=None):
    # min_duration: [min on duration, min off duration]

    # 1. OFF-ON transitions
    _, valid = fma.general.restrict(on[:,0]-.0001,isa,s_ind=True)
    if min_duration is not None:
        valid_where = np.where(valid)[0]
        valid_for_off = (on[valid_where,0] - on[valid_where-1,1] > min_duration[1]) # duration of preceding OFF
        valid_for_on = (on[valid_where,1] - on[valid_where,0] > min_duration[0]) # duration of following ON
        valid = valid_where[valid_for_off & valid_for_on]
    off_on = on[valid,0]

    # 1. ON-OFF transitions
    _, valid = fma.general.restrict(on[:,1]+.0001,isa,s_ind=True)
    if min_duration is not None:
        print('implement!') # WATCH OUT FOR ASSUMPTION THAT isa ENDS WITH ON
    on_off = on[valid,1]

    return off_on, on_off


def ISPhase(phase_times,dt=0.0001,isa=None):

    range = np.array([0,2*np.pi])
    t_phi = np.concatenate(list(phase_times.values()))
    phi = np.concatenate([np.full(len(times),phase) for phase, times in phase_times.items()])

    # sort by time
    t_phi, valid = np.unique(t_phi,return_index=True)
    phi = phi[valid]

    # unwrap phase
    unwrap_idx = phi[1:] < phi[:-1]
    unwrap_idx = np.insert(np.cumsum(unwrap_idx),0,0)
    phi += unwrap_idx * np.diff(range)

    # interpolate and wrap phase
    t = np.arange(t_phi[0], t_phi[-1]+dt/2, dt)
    phi = np.interp(t,t_phi,phi)
    phi = np.mod(phi,np.diff(range))
    phi = np.column_stack((t,phi))

    # restrict to ISA
    if isa is not None:
        phi = fma.general.restrict(phi,isa)

    # fraction of time spent in ON
    on_fraction = np.sum((phi[:,1] > range[0]) & (phi[:,1] < range.mean())) / len(phi)

    return phi, on_fraction


def isCoupled(times_a,times_b,windows):
    # smallest times_b following each event a
    idx = np.searchsorted(times_b,times_a,side='right')
    valid = idx < len(times_b)
    distance = np.full_like(times_a,np.nan,dtype=float)
    distance[valid] = times_b[idx[valid]] - times_a[valid]
    is_coupled_a = (distance >= windows[0]) & (distance <= windows[1])
    # highest times_a preceding each event b
    idx = np.searchsorted(times_a,times_b,side="left") - 1
    valid = idx >= 0
    distance = np.full_like(times_b,np.nan,dtype=float)
    distance[valid] = times_b[valid] - times_a[idx[valid]]
    is_coupled_b = (distance >= windows[0]) & (distance <= windows[1])
    return is_coupled_a, is_coupled_b


def loadHpcPfcEvents(session,names=None,coupl=None,delta='all',isa=None):

    # 1. load
    if names is None: names = ['ripples','deltaWaves','spindles']
    events = {}
    for name in names:
        try:
            ev, _ = fma.data.loadEvents(session,name)
            events[name] = ev[name]['col1'] if 'col1' in ev[name] else ev[name]['peaks']
        except FileNotFoundError:
            events[name] = []
    out = (events,)

    # 2. assess coupling
    if coupl:
        is_coupled = {}
        is_coupled['ripples'], is_coupled['deltas_rip'] = isCoupled(events['ripples'],events['deltaWaves'],[0.05,0.25]) # [0.05,0.25]
        is_coupled['deltas_spin'], is_coupled['spindles'] = isCoupled(events['deltaWaves'],events['spindles'],[0.1,1.3]) # [0.1,1.3]
        match delta:
            case 'ripples':
                is_coupled['deltaWaves'] = is_coupled['deltas_rip']
            case 'spindles':
                is_coupled['deltaWaves'] = is_coupled['deltas_spin']
            case _:
                is_coupled['deltaWaves'] = is_coupled['deltas_rip'] & is_coupled['deltas_spin']
        out = (events, is_coupled)

    # 3. assess ISA
    if isa is not None:
        is_isa = {e: fma.general.restrict(events[e],isa,s_ind=True)[1] for e in events}
        out = (*out, is_isa)

    return out if len(out) > 1 else out[0]


def epochIndex(events,epochs,duration=None,nrem=None):

    if duration is not None:
        nrem = fma.general.intersectIntervals((nrem,epochs))
        _, nrem_epoch_idx = fma.general.restrict(nrem[:,0],epochs,i_ind=True)

        # cumulative nREM duration per epoch
        nrem_start = np.full(len(epochs),np.nan)
        u, idx = np.unique(nrem_epoch_idx,return_index=True)
        nrem_cum_dur = np.insert(np.cumsum(np.diff(nrem)),0,0)
        nrem_start[u] = nrem_cum_dur[idx]
        nrem_start[np.isnan(nrem_start)] = nrem_start[np.where(np.isnan(nrem_start))[0] + 1]

    epoch_idx = {}
    for event in events:

        # assign epoch indexes
        events[event], epoch_idx[event] = fma.general.restrict(events[event],epochs,i_ind=True)
        epoch_idx[event] += 1

        # shift to compute duration
        if duration is not None:
            shifted, nrem_idx = fma.general.restrict(events[event],nrem,shift=True,i_ind=True)
            shifted -= nrem_start[nrem_epoch_idx[nrem_idx]] # make shifted times relative to their epoch start
            epoch_idx[event][shifted > duration] *= -1

    return events, epoch_idx


''' old functions

# function to train GLM
def _ISRPopulationGLM(session,n_delays=1):

    data = loadMATData(session,'populGLM','PopulGLM')
    data['y'] = data['y'].ravel()[n_delays:-n_delays]
    n_times = data['y'].shape[0]
    prediction = {}
    score = {}
    X_out = {}

    # divide data in test and train set
    ids = np.repeat(range(5),np.ceil((n_times)/5))[:n_times]

    keys = [k for k in data.keys() if k != 'y']
    for key in keys:

        prediction[key] = np.full(data['y'].shape,np.nan)
        score[key] = np.zeros((5,2*n_delays-1))

        for j in range(-n_delays+1,n_delays):

            # apply delay
            X_delayed = data[key][j+n_delays:j+n_delays+n_times,:]
            # output non-delayed X
            if j == 0:
                X_out[key] = X_delayed
            for i in range(5):

                is_train = ids != i
                mdl = Pipeline([
                    ('interactions', PolynomialFeatures(degree=2,interaction_only=True)),
                    ('scaler', StandardScaler()),
                    ('elastic_net', ElasticNetCV())])
                mdl.fit(X_delayed[is_train,:],data['y'][is_train])
                score[key][i,j+n_delays-1] = mdl.score(X_delayed[~is_train,:],data['y'][~is_train])
                if j == 0:
                    prediction[key][~is_train] = mdl.predict(X_delayed[~is_train,:])

    return data['y'], prediction, score, X_out

# try Gaussian Process Regression
def _ISRPopulationGPR(session):

    data = loadMATData(session,'populGLM','PopulGLM')
    n_delays = 0 #10
    data['y'] = data['y'].ravel()[n_delays:] # WAS :-n_delays
    n_times = data['y'].shape[0]
    prediction = {}
    std_prediction = {}
    score = {}
    temp = {}
    kernel = 1 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e2)) # + WhiteKernel(noise_level=1.0)

    # divide data in test and train set
    ids = np.repeat(range(5),np.ceil((n_times)/5))[:n_times]

    keys = [k for k in data.keys() if k != 'y']
    for key in keys:

        prediction[key] = np.full(data['y'].shape,np.nan)
        std_prediction[key] = np.full(data['y'].shape,np.nan)
        score[key] = np.zeros((5,1)) # WAS np.zeros((5,2*n_delays-1))

        # standardize predictors
        X = zscore(data[key])
        temp[key] = X

        for j in [0]: # WAS range(-n_delays+1,n_delays):

            X_delayed = X[j+n_delays:j+n_delays+n_times,:]
            for i in range(5):

                is_train = ids != i
                mdl = GaussianProcessRegressor(kernel=kernel,alpha=0.1,n_restarts_optimizer=9)
                mdl.fit(X_delayed[is_train,:],data['y'][is_train])
                prediction[key][~is_train], std_prediction[key][~is_train] = mdl.predict(X_delayed[~is_train,:],return_std=True)
                score[key][i] = mdl.score(X_delayed[~is_train,:],data['y'][~is_train])

    return data['y'], prediction, score, temp
'''