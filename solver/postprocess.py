"""
reservoir simulation project 1(2022)
2D Multiphase reservoir simulation: Postprocessing file
Author: Promise O. Longe
Email: longepromise@ku.edu
Date modified: 04/18/2022
"""

import numpy as np
import matplotlib.pyplot as plt

def postprocess(P_plot,numerical,well,time):
    nmax = int(numerical.tfinal / numerical.dt)
    Pwf  = np.zeros((nmax+1,len(well.x)))
    BlockPress = np.zeros((nmax+1,len(well.x)))
    qwf_w  = np.zeros((nmax+1,len(well.x)))
    qwf_o  = np.zeros((nmax+1,len(well.x)))


    for j in range(0,len(well.x)):
        for i in range(0,nmax):
            if well.type[j][0] == 1:  # 1 for rate
                if well.constraint[j][0] > 0.0:  # injector well
                    well.fw[j][i] = 1.0
                Pwf[i,j]  = well.Pwf[j][i]
                BlockPress[i, j] = well.BlockPress[j][i]
                qwf_w[i,j]  = well.qwf_w[j][i]
                qwf_o[i, j] = 0.0

            elif (well.typetime[i][j] == 2):        # 2 for BHP
                Pwf[i,j]  = well.Pwf[j][i]
                BlockPress[i, j] = well.BlockPress[j][i]
                qwf_o[i,j]  = well.qwf_o[j][i]
                qwf_w[i,j]  = well.qwf_w[j][i]
    print(Pwf)
    print(qwf_o)
    print(qwf_w)

    np.savetxt('blockPress.csv', BlockPress, delimiter=',', fmt=['%f', '%f', '%f'], header='Well(3,2), Well(4,4), Well(9,3)')
    np.savetxt('Pwf.csv', Pwf, delimiter=',', fmt=['%f', '%f', '%f'], header='Well(3,2), Well(4,4), Well(9,3)')
    np.savetxt('qwf_o.csv', qwf_o, delimiter=',', fmt=['%f', '%f', '%f'], header='Well(3,2), Well(4,4), Well(9,3)')
    np.savetxt('qwf_w.csv', qwf_w, delimiter=',', fmt=['%f', '%f', '%f'], header='Well(3,2), Well(4,4), Well(9,3)')

    #print('BHP at 10 days', Pwf[10,:])
    #print('BHP at 20 days', Pwf[20,:])

    #print('Flow rate at 10 days', qwf_o[1,:])
    #print('Flow rate at 20 days', qwf_o[10,:])
    
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plt.plot(time[0:len(time)-1],(qwf_w[0:len(time)-1,0]),'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],(qwf_w[0:len(time)-1,1]),'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],(qwf_w[0:len(time)-1,2]),'b',label=f'Well 3')
    manager = plt.get_current_fig_manager()
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Water Rate $[bbl/day]$')
    plt.xlabel(r'time [days]')
    plt.savefig('WaterRatevstime.png',bbox_inches='tight')

    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plt.plot(time[0:len(time)-1],(qwf_o[0:len(time)-1,0]),'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],(qwf_o[0:len(time)-1,1]),'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],(qwf_o[0:len(time)-1,2]),'b',label=f'Well 3')
    manager = plt.get_current_fig_manager()
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Oil Rate $[bbl/day]$')
    plt.xlabel(r'time [days]')
    plt.savefig('OilRatesvstime.png',bbox_inches='tight')
    
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,0],'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,1],'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],Pwf[0:len(time)-1,2],'b',label=f'Well 3')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    manager = plt.get_current_fig_manager()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'BHP [psi]')
    plt.xlabel(r'time [days]')
    plt.savefig('BHPvstime.png',bbox_inches='tight')

    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plt.plot(time[0:len(time)-1],BlockPress[0:len(time)-1,0],'k',label=f'Well 1')
    plt.plot(time[0:len(time)-1],BlockPress[0:len(time)-1,1],'r',label=f'Well 2')
    plt.plot(time[0:len(time)-1],BlockPress[0:len(time)-1,2],'b',label=f'Well 3')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    manager = plt.get_current_fig_manager()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Grid block Pressure [psi]')
    plt.xlabel(r'time [days]')
    plt.savefig('BlockPress.png',bbox_inches='tight')
    
    fig = plt.figure(figsize=(15,7.5) , dpi=100)
    plt.plot(time[0:len(time)-1],well.fw[0,0:len(time)-1],'k-',label=f'Well 1')
    plt.plot(time[0:len(time)-1],well.fw[1,0:len(time)-1],'r-',label=f'Well 2')
    plt.plot(time[0:len(time)-1],well.fw[2,0:len(time)-1],'b-',label=f'Well 3')
    plt.legend(loc='best', shadow=False, fontsize='x-large')
    manager = plt.get_current_fig_manager()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.ylabel(r'Water cut $f_w$')
    plt.xlabel(r'time [days]')
    plt.savefig('Watercutvstime.png',bbox_inches='tight')
    
    return