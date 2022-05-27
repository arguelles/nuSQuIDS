import proposal as pp
import numpy as np
import h5py

medium = pp.medium.StandardRock()
energy_cuts = pp.EnergyCutSettings(-1, -1)

print(medium)

interpolation_def = pp.InterpolationDef()
interpolation_def.path_to_tables = "/n/home05/agarciasoto/software/proposal/outputs/proposal_tables"

cut = 0.001

particle_def = pp.particle.TauMinusDef()
brems = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
pair  = pp.parametrization.pairproduction.KelnerKokoulinPetrukhin(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
photo = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(shadow_effect=pp.parametrization.photonuclear.ShadowDuttaRenoSarcevicSeckel(),particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
ion   = pp.parametrization.ionization.BetheBlochRossi(particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)

file_h5 = h5py.File('./tau_rock.h5', 'w')

energies = np.logspace(2,11,91)
ys       = np.logspace(-20,0,10001)

xsec    = np.zeros(len(energies))
beta    = np.zeros(len(energies))
betacut = np.zeros(len(energies))
dsdy    = np.zeros((len(ys),len(energies)))

for i in range(len(energies)) :            
    for j in range(len(ys[:-1])) :
        aux_brems = brems.differential_crosssection(energies[i]*1e3,ys[j])
        aux_pair  = pair.differential_crosssection(energies[i]*1e3,ys[j])
        aux_photo = photo.differential_crosssection(energies[i]*1e3,ys[j])
        aux_ion   = ion.differential_crosssection(energies[i]*1e3,ys[j])
        if np.isnan(aux_brems) or aux_brems<0 : aux_brems = 0
        if np.isnan(aux_pair) or aux_pair<0  : aux_pair  = 0
        if np.isnan(aux_photo) or aux_photo<0 : aux_photo = 0
        if np.isnan(aux_ion) or aux_ion<0   : aux_ion   = 0            
        aux = aux_brems + aux_pair + aux_photo + aux_ion

        if aux==0 : aux = 1e-50

        dsdy[j][i] = aux
        width = (ys[j+1]-ys[j]) 
        xsec[i] += dsdy[j][i]*width
        beta[i] += ys[j]*dsdy[j][i]*width   
        if ys[j]<0.001 : betacut[i] += ys[j]*dsdy[j][i]*width   

    aux_brems = brems.differential_crosssection(energies[i]*1e3,ys[-1])
    aux_pair  = pair.differential_crosssection(energies[i]*1e3,ys[-1])
    aux_photo = photo.differential_crosssection(energies[i]*1e3,ys[-1])
    aux_ion   = ion.differential_crosssection(energies[i]*1e3,ys[-1])
    if np.isnan(aux_brems) or aux_brems<0 : aux_brems = 0
    if np.isnan(aux_pair) or aux_pair<0  : aux_pair  = 0
    if np.isnan(aux_photo) or aux_photo<0 : aux_photo = 0
    if np.isnan(aux_ion) or aux_ion<0   : aux_ion   = 0            
    aux = aux_brems + aux_pair + aux_photo + aux_ion
    if aux==0 : aux = 1e-50
    dsdy[-1][i] = aux

energies = np.log10(energies*1e9)
ys       = ys
xsec     = np.log10(xsec/medium.sum_nucleons)
beta     = np.log10(beta/medium.sum_nucleons)
betacut  = np.log10(betacut/medium.sum_nucleons)
dsdy     = np.log10(dsdy/medium.sum_nucleons)

file_h5.create_dataset( "energies", data=energies )
file_h5.create_dataset( "ys",       data=ys       )
file_h5.create_dataset( "s",        data=xsec     )
file_h5.create_dataset( "beta",     data=beta     )
file_h5.create_dataset( "betacut",  data=betacut  )
file_h5.create_dataset( "dsdy",     data=dsdy     )

file_h5.close()


'''
ion
    models['BetheBlochRossi']          = pp.parametrization.ionization.BetheBlochRossi(particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['IonizBergerSeltzerBhabha'] = pp.parametrization.ionization.BergerSeltzerBhabha(particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['IonizBergerSeltzerMoller'] = pp.parametrization.ionization.BergerSeltzerMoller(particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)

pair
    models['KelnerKokoulinPetrukhin1']   = pp.parametrization.pairproduction.KelnerKokoulinPetrukhin(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['KelnerKokoulinPetrukhin2']   = pp.parametrization.pairproduction.KelnerKokoulinPetrukhin(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['SandrockSoedingreksoRhode1'] = pp.parametrization.pairproduction.SandrockSoedingreksoRhode(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['SandrockSoedingreksoRhode2'] = pp.parametrization.pairproduction.SandrockSoedingreksoRhode(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)

brems
    models['PetrukhinShestakov1']        = pp.parametrization.bremsstrahlung.PetrukhinShestakov(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['PetrukhinShestakov2']        = pp.parametrization.bremsstrahlung.PetrukhinShestakov(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['KelnerKokoulinPetrukhin1']   = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['KelnerKokoulinPetrukhin2']   = pp.parametrization.bremsstrahlung.KelnerKokoulinPetrukhin(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['CompleteScreening1']         = pp.parametrization.bremsstrahlung.CompleteScreening(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['CompleteScreening2']         = pp.parametrization.bremsstrahlung.CompleteScreening(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['AndreevBezrukovBugaev1']     = pp.parametrization.bremsstrahlung.AndreevBezrukovBugaev(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['AndreevBezrukovBugaev2']     = pp.parametrization.bremsstrahlung.AndreevBezrukovBugaev(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['SandrockSoedingreksoRhode1'] = pp.parametrization.bremsstrahlung.SandrockSoedingreksoRhode(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['SandrockSoedingreksoRhode2'] = pp.parametrization.bremsstrahlung.SandrockSoedingreksoRhode(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['ElectronScreening1']         = pp.parametrization.bremsstrahlung.ElectronScreening(lpm_effect=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['ElectronScreening2']         = pp.parametrization.bremsstrahlung.ElectronScreening(lpm_effect=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)

photo
    models['Zeus1']                      = pp.parametrization.photonuclear.Zeus(add_pertubative=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1.)
    models['Zeus2']                      = pp.parametrization.photonuclear.Zeus(add_pertubative=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1.)
    models['BezrukovBugaev1']            = pp.parametrization.photonuclear.BezrukovBugaev(add_pertubative=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1.)
    models['BezrukovBugaev2']            = pp.parametrization.photonuclear.BezrukovBugaev(add_pertubative=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1.)
    models['Rhode1']                     = pp.parametrization.photonuclear.Rhode(add_pertubative=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['Rhode2']                     = pp.parametrization.photonuclear.Rhode(add_pertubative=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['Kokoulin1']                  = pp.parametrization.photonuclear.Kokoulin(add_pertubative=False,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['Kokoulin2']                  = pp.parametrization.photonuclear.Kokoulin(add_pertubative=True,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['AbramowiczLevinLevyMaor911'] = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor91(shadow_effect=shadow1,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['AbramowiczLevinLevyMaor912'] = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor91(shadow_effect=shadow2,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['AbramowiczLevinLevyMaor971'] = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(shadow_effect=shadow1,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['AbramowiczLevinLevyMaor972'] = pp.parametrization.photonuclear.AbramowiczLevinLevyMaor97(shadow_effect=shadow2,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['ButkevichMikhailov1']        = pp.parametrization.photonuclear.ButkevichMikhailov(shadow_effect=shadow1,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['ButkevichMikhailov2']        = pp.parametrization.photonuclear.ButkevichMikhailov(shadow_effect=shadow2,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['RenoSarcevicSu1']            = pp.parametrization.photonuclear.RenoSarcevicSu(shadow_effect=shadow1,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
    models['RenoSarcevicSu2']            = pp.parametrization.photonuclear.RenoSarcevicSu(shadow_effect=shadow2,particle_def=particle_def,medium=medium,energy_cuts=energy_cuts,multiplier=1)
'''
