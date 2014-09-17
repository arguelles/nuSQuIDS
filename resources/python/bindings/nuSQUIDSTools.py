import numpy as np
import matplotlib.pyplot as plt
import nuSQUIDSpy as nsq

# we will use the boost injector metaclass

BoostPythonMetaclass = nsq.nuSQUIDS.__class__

class injector(object):
    class __metaclass__(BoostPythonMetaclass):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)
            return type.__init__(self, name, bases, dict)

# nuSQUIDS python extensions

class ExtNuSQUIDS(injector,nsq.nuSQUIDS):
    """ nuSQUIDS class extension """

    def __SUBasis(i,nsun):
        return 0


    def PrintStatus(self,filepath):
        """ Print the nuSQUIDS status to a file"""
        file = open(filepath,"w")
        file.write('% -- nuSQuIDS Automatically Generated Report -- \n')
        file.write('\\documentclass{article}\n')
        file.write('\\usepackage{fullpage}\n')

        file.write('\\begin{document}\n')

        file.write('\\section{Basic}\n')
        file.write('Number of neutrinos:\t{0}'.format(self.GetNumNeu()))
        file.write('\\section{Physics Parameters}\n')
        file.write('\\subsection{Mixing Angles}\n')
        file.write('\\subsection{CP Phases}\n')
        file.write('\\subsection{Square Mass Differences}\n')
        file.write('\\subsection{Mixing Matrix}\n')
        file.write('\\section{Hamiltonian}\n')

        file.write('\\end{document}\n')
        file.close()

    def PrintStatus(self):
        return 0

    def PlotFlavor(self):
        return 0

    def PlotMass(self):
        return 0

    def PlotComposition(self):
        return 0

    def PrintH0(self):
        return 0

    def PrintHI(self):
        return 0

    def PrintState(self,ie,rho):
        return 0


