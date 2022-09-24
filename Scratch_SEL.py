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

Materials=['Berea_Sandstone_Specimen_2']
markerstyle=['s','o','d','D','p','h','<','>','v','^','1','2','s','o','d']
Colors=['red','blue','green','purple','orangered','green','orange','turquoise','orangered','olive','navy','plum']
Edgecolors=['darkred','darkblue','darkgreen','black','black']




##tests=[filetest for filetest in glob.glob('*.txt') if (material+'#') in filetest]






filetest=['Tests']
filetest_number=[11]  ## number of tests

for j in range(len(filetest)):

    Kcall=[]
    D0all=[]
    xmin,xmax,ymin=1e-1,0,0.1

    fig=plt.figure(1)
    ax = plt.subplot()
    axis_font = {'fontname':'Arial', 'size':'28'}
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontname('Arial')
        label.set_fontsize(24)
    ax.set_position([0.15,0.15,0.8,0.8])
    ax.set_axisbelow(True)

    for n in range(1,filetest_number[j]):


        ####------- colors
        colors_ = list(six.iteritems(colors.cnames))
        # Add the single letter colors.
        for name, rgb in six.iteritems(colors.ColorConverter.colors):
            hex_ = colors.rgb2hex(rgb)
            colors_.append((name, hex_))
        filename=filetest[j]+' '+str(n)+'.TXT'
        filename2=os.path.splitext(filename)[0]
        file=open(filename,"r")
        lines=file.readlines()
        file.close()

        Ft=[float(lines[k+36].split()[3]) for k in range(len(lines)-36)]
        Pd=[float(lines[k+36].split()[4]) for k in range(len(lines)-36)]
        X=[float(lines[k+36].split()[0]) for k in range(len(lines)-36)]
        Pf=[float(lines[k+36].split()[6]) for k in range(len(lines)-36)]

        Ft=numpy.array(Ft)
        Pd=numpy.array(Pd)
        X=numpy.array(X)
        Pf=numpy.array(Pf)
        Ft=Ft-Ft[0]
        Pd=Pd-Pd[0]

        MN=40

        Ft= numpy.convolve(Ft, numpy.ones((MN,))/MN, mode='valid')
        Pd= numpy.convolve(Pd, numpy.ones((MN,))/MN, mode='valid')




        ######################  Spherical Scaling ############################

        find=numpy.where((Pd>0.) & (Pd<27) & (Ft>0)) ## spherical
        Ft1=numpy.array([Ft[k] for k in range(len(Ft)) if (Pd[k]>0) and (Pd[k]<27) and (Ft[k]>0.)])
        Pd1=numpy.array([Pd[k] for k in range(len(Ft)) if (Pd[k]>0) and (Pd[k]<27) and (Ft[k]>0.)])


           
        ##---- calculate SEl parameters
        R=200 ## tip radius in nm
        Reff=54.51*R/(32.0/3.0) ### effective tip radius

        A=1e-6*(Reff**2)*2/3*((2*Pd1/Reff)**(1.5)) ### area in mm^2

        sigma_N=Ft1/(1000*A)  ### nominal stress in GPa
        xdata=numpy.log(Pd1/3)
        ydata=sigma_N


        ## curve fitting to logarithmic creep
        def fcn2min(params, x, data):
            Kc = params['Kc']
            D0 = params['D0']
            
            return ydata-Kc/numpy.sqrt(D0+numpy.exp(x))
            
        params = Parameters()
        params.add('Kc', value=0.8, min=0.1,max=1.8)
        params.add('D0',min=0.1,value=2, max=9)


        minner = Minimizer(fcn2min, params, fcn_args=(xdata, ydata))
        result1 = minner.minimize(method)

        Kc=result1.params.valuesdict()['Kc']  
        D0=result1.params.valuesdict()['D0']  ## D0 in micrometer
        Bft=Kc/numpy.sqrt(D0)

        fi=Bft*((1+numpy.exp(xdata)/(D0))**(-0.5))
        SStot=sum((sigma_N-numpy.mean(sigma_N))**2)
        SSres=sum((sigma_N-fi)**2)
        ##R2=1 - result1.residual.var() / numpy.var(ydata)
        R2=1-SSres/SStot
        ##ci = lmfit.conf_interval(minner, result) ### confidence intervals
        RMSE=numpy.mean(result1.residual*result1.residual)


        ####------------------ plot


        print( filename, R2, RMSE, Kc)

        Kcall.append(Kc)
        D0all.append(D0)

        plt.plot(Pd1/(3*D0),sigma_N/Bft, '--'+markerstyle[n],color=str(n*1.0/(filetest_number[j]+2)),markeredgecolor='black', ms=8, markeredgewidth=0.05)
        xmax=max(xmax,max(Pd1/(3*D0)))
        xmin=min(xmin, min(Pd1/(3*D0)))
        ymin=min(ymin,min(sigma_N/Bft))



        ##t=plt.text(xmin*1.01, 0.3+0.07,'$K_c$='+"%10.2f" %numpy.mean(numpy.array(Kcall))+'$\pm$'+"%10.2f" %numpy.std(numpy.array(Kcall)) +' MPa.m$^{0.5}$',fontsize=24)
        ##t.set_bbox(dict(facecolor='white',edgecolor='white'))
        ##plt.text(xmin*1.01, 0.3+0.03,'$D_0$='+"%10.2f" %numpy.mean(numpy.array(D0all))+'$\pm$'+"%10.2f" %numpy.std(numpy.array(D0all))+ ' $\mu$m',fontsize=24)
        ##t.set_bbox(dict(facecolor='white', edgecolor='white'))

    aux=numpy.array([10.0**i for i in numpy.linspace(-4,2,100)])

    ytheo=(1+aux)**(-0.5)
    plt.plot(aux,ytheo,'-k',linewidth=2)
    ax.set_ylim(0.,2)
    ax.set_xlim(0,100)
    ax.set_yticks([0.3,0.5,0.7,0.9,1.,1.2,1.4,1.6])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlabel("$D/D0$", **axis_font)
    plt.ylabel('$\sigma_N/Bf^\prime_t$', **axis_font)
    plt.grid(True,which="both",color='lightgrey')
    plt.savefig('Toughness_'+filetest[j]+'_SEL_v2' + '.png')

    plt.show()
    plt.close(fig)


    ##
    ##



