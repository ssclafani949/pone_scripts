import tables
import numpy as np
import histlite as hl
import pandas as pd
import matplotlib.pyplot as plt

leptoninj_high_stats = tables.open_file('../data/LI_ranged_high_stats.hdf5')



def get_weights(f, nugen=True, nfiles=100):
    #print('nugen', nugen)
    #print('nfiles : {}'.format(nfiles) )
    event = f.root.I3EventHeader.cols.Event[:]
    nevents = 100
    #print(len(event)/nfiles)
    #print(nevents)
    if nugen:
        true_zenith = f.root.NuGPrimary.cols.zenith[:]
        true_azimuth = f.root.NuGPrimary.cols.azimuth[:]
        true_energy = f.root.NuGPrimary.cols.energy[:]
        #print(min(true_energy))
        oneweight = f.root.I3MCWeightDict.cols.OneWeight[:]
        ow = oneweight / (nfiles * nevents)
        w_astro = 1e-18 * ow * (true_energy/1e5)**-2
        array = np.vstack([ true_zenith, true_azimuth, true_energy,  ow, w_astro])
        df = pd.DataFrame(array.T, 
            columns= ['true_zenith', 'true_azimuth', 'true_energy',  'ow', 'w_astro'])
    else:
        true_zenith = f.root.EventProperties.cols.zenith[:]
        true_azimuth = f.root.EventProperties.cols.azimuth[:]
        true_energy = f.root.EventProperties.cols.totalEnergy[:]
        flux = f.root.LeptonInjection_flux.cols.value[:]
        weight = f.root.LeptonInjection_weight.cols.value[:]
        oneweight = f.root.LeptonInjection_oneweight.cols.value[:]
        surv_prob = f.root.LeptonInjection_SurvivalProb.cols.value[:]
        ow = (oneweight * surv_prob) / (nfiles)
        w_astro = 1e-18 * ow * (true_energy/1e5)**-2
        array = np.vstack([ true_zenith, true_azimuth, true_energy, weight, ow, w_astro])
        df = pd.DataFrame(array.T, 
            columns= ['true_zenith', 'true_azimuth', 'true_energy',  'weight', 'ow', 'w_astro'])
    return df




def effa_plot_notlogged(ax, sim, key, zen_deg_min = 0, zen_deg_max = 180, 
                                selection_eff = 1, alpha=1, label = 'MISC', color = 'k', ls = '-'):
    a = sim[key]
    mask = (a.true_zenith < np.radians(zen_deg_max)) & (a.true_zenith > np.radians(zen_deg_min))

    dlogE=.1
    bins = 10**np.arange(3,7.01,dlogE)
    solid_angle= 2*np.pi* (-1*(np.cos(np.radians(zen_deg_max))) - -1*np.cos(np.radians(zen_deg_min)))
    print(np.degrees(solid_angle))
    area=  (1/(1e4*np.log(10)) * (a.ow[mask] / ( (a.true_energy[mask]) *  solid_angle *  dlogE)))
    h = hl.hist((a.true_energy[mask]), 
        weights = selection_eff * area, bins=bins);
    hl.plot1d(ax, h, histtype='step', linewidth=2, color=color, alpha=alpha, label=label, ls = ls, log=True)
    ax.loglog()
    ax.grid()
    ax.set_ylabel('$A_\mathsf{Eff}$($m^2$)')
    ax.set_xlabel('log$_{10}$(E) [1/GeV]')
    plt.tight_layout()
    return h



sim= {}
names = ['leptoninj_highstats']

nugen = [False]
nfiles = [1000]


for i, (name, mc) in enumerate(zip(names, [leptoninj_high_stats])):
    sim[name] = get_weights(mc, nugen=nugen[i], nfiles=nfiles[i])
    mask = (sim[name].true_energy > 1000) & (sim[name].true_energy < 1e7)
    print('{} astro events / year: {:.3}'.format(name, sum(sim[name][mask].w_astro) * 86400 * 365))


fig,(ax1) = plt.subplots(1,1, figsize=(8,6))

zen_bands = np.arange(0,181,30)
colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
h_all = effa_plot_notlogged(ax1, sim,  key='leptoninj_highstats', label='Allsky', color='blue', alpha=1, ls='--')
for i, zen_min in enumerate(zen_bands):
    if zen_min < 180:
        print(zen_min)
        h_l = effa_plot_notlogged(ax1, sim,  zen_deg_min=zen_bands[i], zen_deg_max = zen_bands[i+1], 
                key='leptoninj_highstats', label='{:.2f} < zen <{:.2f}'.format(zen_bands[i], zen_bands[i+1]), 
                color=colors[i], alpha=1, ls='-')

icecube_effa = [55.059377358189856, 0.0008646535029500378,
        83.9857063179509, 0.0028700848540715793,
        179.58858500014657, 0.01767536622987677,
        406.25802790141057, 0.08751958772041732,
        868.7109696923557, 0.37470044642998684,
        2199.369957682293, 1.6636142493842194,
        7378.54380841416, 8.237387069570994,
        16691.4431470578, 21.983926488622846,
        42258.77177653403, 52.608101759655035,
        104019.53653649216, 94.12049672680651,
        424973.5083400382, 194.74830399087512,
        1347672.455433188, 348.42175435383876,
        4039773.0556262657, 449.39845907216517,
        16046450.217745425, 720.93272022235]
icecube_energy = icecube_effa[::2]
icecube_effa_m2 = icecube_effa[1::2]

icecube_effa2 = [83.99, 0.01,
        151.68, 0.04,
        333.60, 0.13,
        1088.12, 1.20,
        6231.91, 12.29,
        17658.09, 37.93,
        54442.91, 97.61,
        303153.50, 312.42,
        1425719.96, 775.31,
        5353112.33, 1438.45,
        14747016.17, 2767.61,
        17460383.96, 2767.61]
icecube2_energy = icecube_effa2[::2]
icecube2_effa_m2 = icecube_effa2[1::2]

#ax1.plot(icecube_energy, icecube_effa_m2, label='IceCube Trigger - IC86', c='brown')
#ax1.plot(icecube2_energy, icecube2_effa_m2, label='IceCube Trigger - AllSky numu + nutau', c='green', ls='-')

ax1.legend(ncol=2)

plt.savefig('/data/p-one/ssclafani/effa_bands.png')

