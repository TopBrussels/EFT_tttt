import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM

def main():
    np.set_printoptions(edgeitems=3)
    np.core.arrayprint._line_width = 120

    eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)

    #Plots with limit contours
    make_plot1d(eft.S)
    #make_plot2d(eft.S)

def C0(C0,s):
        '''
        Calculations of the tttt EFT cross section based on just one operator C0
        :param C0:
        :param s:
        :return: sigma_tttt = sigma_tttt_SM +  C_0*sigma_0 + C_0*C_0*sigma_00
        '''
        c = [C0,0.,0.,0.,0.]
        row = [ sig_SM, c[0]*s[0], c[1]*s[1], c[2]*s[2], c[3]*s[3], c[4]*s[4], 
                (c[0]**2.)*s[5], (c[1]**2.)*s[6], (c[2]**2.)*s[7], (c[3]**2.)*s[8], (c[4]**2.)*s[9], 
                2.*c[0]*c[1]*s[10], 2.*c[0]*c[2]*s[11], 2.*c[0]*c[3]*s[12], 2.*c[0]*c[4]*s[13], 
                2.*c[1]*c[2]*s[14], 2.*c[1]*c[3]*s[15], 2.*c[1]*c[4]*s[16], 
                2.*c[2]*c[3]*s[17], 2.*c[2]*c[4]*s[17], 
                2.*c[3]*c[4]*s[19]]
        return sum(row)

def C0C1(C0,C1,s):
        '''
        Calculations of the tttt EFT cross section based on just two operators C0, C1
        :param C0:
        :param s:
        :return: sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij, where i,j = 0,1
        '''
        c = [C0,C1,0.,0.,0.]
        row = [ sig_SM, c[0]*s[0], c[1]*s[1], c[2]*s[2], c[3]*s[3], c[4]*s[4], 
                (c[0]**2.)*s[5], (c[1]**2.)*s[6], (c[2]**2.)*s[7], (c[3]**2.)*s[8], (c[4]**2.)*s[9], 
                2.*c[0]*c[1]*s[10], 2.*c[0]*c[2]*s[11], 2.*c[0]*c[3]*s[12], 2.*c[0]*c[4]*s[13], 
                2.*c[1]*c[2]*s[14], 2.*c[1]*c[3]*s[15], 2.*c[1]*c[4]*s[16], 
                2.*c[2]*c[3]*s[17], 2.*c[2]*c[4]*s[17], 
                2.*c[3]*c[4]*s[19]]
        return sum(row)

def make_plot1d(sigma):
        xlist = np.linspace(-3.0, 3.0, num=120)
        f = np.vectorize( C0, excluded=set([1]))
        ylist = f(xlist,sigma)

        plt.rc('text', usetex=True)
        plt.rcParams.update({'font.size': 20})
        plt.figure()
        plt.title('')
        plt.plot(xlist,ylist,linewidth=2.0)

        xband = np.linspace(-3.0, 3.0, 2)
        yband = np.array([2.7*sig_SM, 2.7*sig_SM])
        band_half_w_1s_pos = np.array([(1.6)*sig_SM,(1.6)*sig_SM])
        band_half_w_1s_neg = np.array([(1.)*sig_SM,(1.)*sig_SM])
        band_half_w_2s_pos = np.array([(3.9)*sig_SM,(3.9)*sig_SM])
        band_half_w_2s_neg = np.array([(1.52)*sig_SM,(1.52)*sig_SM])

        plt.plot(xband, yband, 'k')
        plt.fill_between(xband, yband-band_half_w_2s_neg, yband+band_half_w_2s_pos, facecolor='green')
        plt.fill_between(xband, yband-band_half_w_1s_neg, yband+band_half_w_1s_pos, facecolor='yellow')

        plt.xlabel(r'$c_{O_{R}}$', labelpad=2)
        plt.ylabel(r'$\sigma_{t\bar{t}t\bar{t}}$ (fb)')

        two_s_band = mpatches.Patch(color='green', label=r'2 s.d.')
        one_s_band = mpatches.Patch(color='yellow', label=r'1 s.d.')
        limit_line = mlines.Line2D([], [], color='black', label='Combined')
        eft_line   = mlines.Line2D([], [], color='blue', label='EFT')
        plt.legend([eft_line,limit_line, one_s_band, two_s_band],['EFT','Combined','1 s.d.','2 s.d.'],loc=4,prop={'size':13})

        plt.show()

def make_plot2d(sigma):
        xlist = np.linspace(-3.0, 3.0, num=120)
        ylist = np.linspace(-3.0, 3.0, num=120)
        X, Y = np.meshgrid(xlist, ylist)
        f = np.vectorize( C0C1, excluded=set([2]))
        Z = f(X,Y,sigma)

        plt.rc('text', usetex=True)
        plt.rcParams.update({'font.size': 20})
        plt.figure()
        levels = [10.0, 20., 30., 40., 50., 70.]
        contour = plt.contour(X, Y, Z, levels, colors='k')
        plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
        level_2015 = [6.4*sig_SM] 
        level_2016_sl = [2.7*sig_SM] 
        contour_2015 = plt.contour(X, Y, Z, level_2015, colors='r',linewidths=np.arange(3.9, 4, .5),linestyles='dashed')
        plt.clabel(contour_2015, colors = 'r', fmt = '2015', fontsize=12)
        contour_2016sl = plt.contour(X, Y, Z, level_2016_sl, colors='r',linewidths=np.arange(3.9, 4, .5))
        plt.clabel(contour_2016sl, colors = 'r', fmt = '2016', fontsize=12)
        contour_filled = plt.contourf(X, Y, Z, 100)
        plt.colorbar(contour_filled)
        plt.title('')
        plt.xlabel(r'$O_{R}$')
        plt.ylabel(r'$O_{L}^{(1)}$')
        plt.show()
    
if __name__ == "__main__":
        main()
