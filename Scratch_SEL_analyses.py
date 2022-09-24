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
When=['Before', 'After',]
When=['After']
Materials=['13-10','13-14']
##Materials=['13-14']
Directions=['parallel','perpendicular']
##Directions=['perpendicular']
markerstyle=['s','o','d','D','p','h','<','>','v','^','1','2']
Colors=['red','blue','green']
Edgecolors=['darkred','darkblue','darkgreen']


file1=open('Scratch_Tests_SEL.out','w')
file1.write('When \t Material \t Direction \t Kc (MPa sqrt(m)) \t D0 (micron) \t R2 \n')



m=0
for material in Materials:
    for when in When:
        for direction in Directions:
            print(when,material, direction)

            tests=[filetest for filetest in glob.glob('*.txt') if (when+'_'+material+'_'+direction) in filetest]
            Kcall=[]
            D0all=[]
            xmin,xmax,ymin=1e-1,0,0.1

              

            for j in range(len(tests)):

                filename=tests[j]
                print(j)
                
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
                Pd=[float(line.split()[5])*1e6 for line in lines[1:-3]] ### penetration depth in microns

                Ft=numpy.array(Ft)
                Pd=numpy.array(Pd)
                X=numpy.array(X)
                Fn=numpy.array(Fn)

               


                
                find=numpy.where((Pd>max(Pd)*0.05) &  (Ft>0))
                Ft=Ft[find]
                Pd=Pd[find]

                MN=20

                Ft= numpy.convolve(Ft, numpy.ones((MN,))/MN, mode='valid')
                Pd= numpy.convolve(Pd, numpy.ones((MN,))/MN, mode='valid')



                ##---- calculate SEl parameters
                R=200 ## tip radius in nm
                theta=numpy.pi/3 ### half-apex angle

                A=1e-6*(Pd**2)*numpy.tan(theta) ### area in mm^2

                sigma_N=Ft/(1000*A)  ### nominal stress in GPa
                xdata=numpy.log(0.25*Pd*numpy.sin(theta))
                ydata=numpy.log(sigma_N)

                
                ## curve fitting to logarithmic creep
                def fcn2min(params, x, data):
                    Kc = params['Kc']
                    D0 = params['D0']
                    return ydata-numpy.log(Kc/numpy.sqrt(D0+numpy.exp(x)))

                params = Parameters()
                params.add('Kc', value=0.6, min=0)
                params.add('D0', value=10, min=0.1)

                minner = Minimizer(fcn2min, params, fcn_args=(xdata, ydata))
                result = minner.minimize(method='tnc')

                Kc=result.params.valuesdict()['Kc']   ### Kc in MPa sqrt(m)
                D0=result.params.valuesdict()['D0']  ## D0 in micrometer
                Bfprime=Kc/numpy.sqrt(D0)
                R21=1 - result.residual.var() / numpy.var(ydata)

##
##                fi=Bfprime*((1+numpy.exp(xdata)/D0)**(-0.5))
##                SStot=sum((sigma_N-numpy.mean(sigma_N))**2)
##                SSres=sum((sigma_N-fi)**2)
##                ##R2=1 - result1.residual.var() / numpy.var(ydata)
##                R2=1-SSres/SStot
##                ##ci = lmfit.conf_interval(minner, result) ### confidence intervals
##                RMSE=numpy.sqrt(numpy.mean((sigma_N-fi)**2/sigma_N**2))
##                
##                ##ci = lmfit.conf_interval(minner, result) ### confidence intervals
##
                print(Kc, R21)
                file1.write(when+'\t'+material+'\t'+direction+'\t'+"%10.2f" % (Kc)+'\t'+"%10.2f" % (D0)+'\t'+"%10.4f" % (R21)+'\n')
 


                   
              
             
##
##                plt.figure(m)
##                ax = plt.subplot()
##                axis_font = {'fontname':'Arial', 'size':'24'}
##                for label in (ax.get_xticklabels() + ax.get_yticklabels()):
##                    label.set_fontname('Arial')
##                    label.set_fontsize(20)
##
##
##                plt.loglog(0.25*Pd*numpy.sin(theta)/D0,sigma_N/Bfprime, '--o', color='r', markerfacecolor='r',markeredgecolor='k', ms=8, markeredgewidth=0.05)
##
##                
##                
##                aux=numpy.array([10.0**i for i in numpy.linspace(-2,5,1000)])
##                plt.loglog(aux,(1+aux)**(-0.5),'k',linewidth=2)
##                ax.set_ylim(0.001,  1.6)
##                ax.set_yticks([0.001,0.01,0.05,0.1,0.3,0.5,0.7,0.9,1.,1.5])
##                ax.set_xlim(0.01,10000)
##                plt.xlabel("$D/D0$", **axis_font)
##                plt.ylabel('$\sigma_N/Bf^\prime_t$', **axis_font)
##                plt.grid(True,which="both",color='lightgrey')
##                plt.tight_layout()
##                plt.savefig('SEL_'+when+'_'+material+'_'+direction +'_'+'%d' %(j)+ '.pdf')
##                plt.close()
##                m=m+1

##
##                plt.figure(m)
##                ax = plt.subplot()
##                axis_font = {'fontname':'Arial', 'size':'24'}
##                for label in (ax.get_xticklabels() + ax.get_yticklabels()):
##                    label.set_fontname('Arial')
##                    label.set_fontsize(20)
##
##
##                plt.plot(0.25*Pd*numpy.sin(theta)/D0,sigma_N/Bfprime, marker='o', markerfacecolor='r',markeredgecolor='k', ms=8, markeredgewidth=0.05)
##
##                
##                
##                aux=numpy.array([10.0**i for i in numpy.linspace(-2,5,1000)])
##                xx=0.25*Pd*numpy.sin(theta)/D0
##                plt.plot(aux,(1+aux)**(-0.5),'k',linewidth=2)
##                ax.set_ylim(0.0,  1.)
##                print(min(0.25*Pd*numpy.sin(theta)/D0),max(0.25*Pd*numpy.sin(theta)/D0))
##                ax.set_xticks([30, 100, 500, 800, 1000])
##                ax.set_xlim(0.0,max(xx)*1.1)
##                plt.xlabel("$D/D0$", **axis_font)
##                plt.ylabel('$\sigma_N/Bf^\prime_t$', **axis_font)
##                plt.grid(True,which="both",color='lightgrey')
##                plt.tight_layout()
##                plt.savefig('SEL_Linear_'+when+'_'+material+'_'+direction +'_'+'%d' %(j)+ '.pdf')
##                plt.close()
##                m=m+1
                
                
file1.close()
                
                
                                        

                


               
                
              

