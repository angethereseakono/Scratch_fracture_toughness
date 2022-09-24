import math
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.ticker
import numpy
import six
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, Model
from scipy.optimize import curve_fit
from lmfit.printfuncs import *
import lmfit
import os
import glob, os




method='tnc'
When=['After']
Materials=['13-10']
Directions=['parallel']
markerstyle=['s','o','d','D','p','h','<','>','v','^','1','2']
Colors=['red','blue','green']
Edgecolors=['darkred','darkblue','darkgreen']
##
##file1=open('All_Scratch_Tests.out','w')
##file1.write('When \t Material \t Direction \t Kc (MPa sqrt(m)) \t dKc MPa sqrt(m) \n')


n=0
for when in When:
    for material in Materials:
        for direction in Directions:
  

            tests=[filetest for filetest in glob.glob('*.txt') if (when+'_'+material+'_'+direction) in filetest]
            Kcall=[]
            D0all=[]
            xmin,xmax,ymin=1e-1,0,0.1

         

            for j in range(len(tests)):

                filename=tests[j]
                
                
                ####------- colors
                colors_ = list(six.iteritems(colors.cnames))
                # Add the single letter colors.
                for name, rgb in six.iteritems(colors.ColorConverter.colors):
                    hex_ = colors.rgb2hex(rgb)
                    colors_.append((name, hex_))

                filename2=os.path.splitext(filename)[0]
                file=open(filename,"r")
                lines=file.readlines()
                file.close()

               
                X=[float(line.split()[0]) for line in lines[1:-3]]
                Ft=[float(line.split()[1]) for line in lines[1:-3]] ### horizontal force in N
                Fn=[float(line.split()[4]) for line in lines[1:-3]] ### vertical force in N
                Pd=[float(line.split()[5])*1e6 for line in lines[1:-3]]#### penetration depth in microns

                ### penetration depth in microns

                Ft=numpy.array(Ft)
                Pd=numpy.array(Pd)
                X=numpy.array(X)
                Fn=numpy.array(Fn)
                Ft=Ft-Ft[0]
                Pd=Pd-Pd[0]

                MN=5

                Ft= numpy.convolve(Ft, numpy.ones((MN,))/MN, mode='valid')
                Pd= numpy.convolve(Pd, numpy.ones((MN,))/MN, mode='valid')

                ##---- calculate fracture toughness
                R=200.0 ## probe tip radius in microns
                f1=13.86*(Pd/R)**3+0*(Pd/R)**2+0.*(Pd/R)

                ###### Normalize the squared horizontal force. gr1=Ft^2/R^3
                gr1=Ft**2./(1e-6*R**3)

                ###perform a linear regression of the normalized squared horizontal force
                ###with the probe shape function
                aux=(gr1/f1)**0.5
                Ks=numpy.mean(aux[Pd>max(Pd)*0.5])
                dKs=numpy.std(aux[Pd>max(Pd)*0.5])
##                file1.write(when+'\t'+material+'\t'+direction+'\t'+"%10.2f" % (Ks)+'\t'+"%10.2f" % (dKs)+'\n')
                print(Ks,dKs)
##
                ##
##                fig=plt.figure(n)
##                ax = plt.subplot()
##                axis_font = {'fontname':'Arial', 'size':'24'}
##                for label in (ax.get_xticklabels() + ax.get_yticklabels()):
##                    label.set_fontname('Arial')
##                    label.set_fontsize(24)
##                ax.set_position([0.15,0.15,0.7,0.8])
##                ax.set_axisbelow(True)
##
##                plt.plot(Pd/R, (gr1/f1)**0.5, linestyle='--', marker='o', color='b')
##                plt.axhline(y = Ks, color = 'k', linestyle = '-')
##                plt.axhline(y = Ks+dKs, color = 'k', linestyle = '--')
##                plt.axhline(y = Ks-dKs, color = 'k', linestyle = '--') 
##                ax.set_ylim(0,7)
##                ax.set_xlim(0,max(Pd/R)*1.05)
##                plt.xlabel("$d/R$", **axis_font)
##                plt.ylabel("$F_T/\sqrt{2pA}$, $MPa\sqrt{m}$", **axis_font)
##                plt.tight_layout()
##                plt.savefig(when+'_'+material+'_'+direction+'_'+'%d' %(j+1) + '.tiff')
##                plt.close()
##                n=n+1
                
                
                

                
file1.close()

               
                
              

