#!/usr/bin/python
"""
Simple script to extract the information
out of the HDF5. It needs a working installation
of PyTables.
"""
import tables as tb
import numpy as np
import sys

if len(sys.argv) != 2 :
    print "Usage: GetTables.py filename.hdf5"

print "Reading file: " + sys.argv[1]
h5file = tb.openFile(sys.argv[1]);
filenameroot = sys.argv[1].rsplit( ".", 1 )[ 0 ]

basic_config = h5file.root.basic;
numneu = h5file.root.basic.getAttr("numneu");
NT = h5file.root.basic.getAttr("NT");

# save flavor composition #
energies = h5file.root.energies;

header_flv = {'neutrino':"E [eV] \t nu_e \t nu_\mu \t nu_\\tau"
         ,'antineutrino':"E [eV] \t nubar_e \t nubar_\mu \t nubar_\\tau"
         ,'both':"E [eV] \t nu_e \t nu_\mu \t nu_\tau \t nubar_e \t nubar_\mu \t nubar_\\tau"}

header_mass = {'neutrino':"E [eV] \t nu_1 \t nu_2 \t nu_3"
         ,'antineutrino':"E [eV] \t nubar_1 \t nubar_2 \t nubar_3"
         ,'both':"E [eV] \t nu_1 \t nu_2 \t nu_3 \t nubar_1 \t nubar_2 \t nubar_3"}

flavortable = [ np.hstack((energies[i],flv)) for i,flv in enumerate(h5file.root.flavorcomp) ]
masstable = [ np.hstack((energies[i],flv)) for i,flv in enumerate(h5file.root.flavorcomp) ]

print "Saving Flavor and Mass composition."
try:
    np.savetxt(filenameroot + "_FlavorComposition.dat", flavortable, header = header_flv[NT])
    np.savetxt(filenameroot + "_MassComposition.dat", masstable, header = header_mass[NT])
except KeyError:
    print "Error: NT not correct; cannot generate header."

print "Saving nuSQuIDS configuration."
