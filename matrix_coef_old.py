import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

sig_SM = 0.009308 * 1000. # fb from MG

def gen_row(c):
        '''
        Generate single row of the matrix for sigma_i and sigma_ij
        EFT predictions has the following analytical expression:
            sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        '''
        row = [c[0], c[1], c[2], c[3], c[4], c[0]**2., c[1]**2., c[2]**2., c[3]**2., c[4]**2., 2.*c[0]*c[1], 2.*c[0]*c[2], 2.*c[0]*c[3], 2.*c[0]*c[4], 2.*c[1]*c[2], 2.*c[1]*c[3], 2.*c[1]*c[4], 2.*c[2]*c[3], 2.*c[2]*c[4], 2.*c[3]*c[4]]
        return row

def gen_eft_xs(c, s):
        '''
        Sum individual contributions from different operators

        :param c: vector of Wilson coefficient values, C_i
        :param s: vector of sigma_i, sigma_ij values
        :return:  sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        '''
        row = [ sig_SM, c[0]*s[0], c[1]*s[1], c[2]*s[2], c[3]*s[3], c[4]*s[4], 
                (c[0]**2.)*s[5], (c[1]**2.)*s[6], (c[2]**2.)*s[7], (c[3]**2.)*s[8], (c[4]**2.)*s[9], 
                2.*c[0]*c[1]*s[10], 2.*c[0]*c[2]*s[11], 2.*c[0]*c[3]*s[12], 2.*c[0]*c[4]*s[13], 
                2.*c[1]*c[2]*s[14], 2.*c[1]*c[3]*s[15], 2.*c[1]*c[4]*s[16], 
                2.*c[2]*c[3]*s[17], 2.*c[2]*c[4]*s[17], 
                2.*c[3]*c[4]*s[19]]
        return sum(row)

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


def main():
    np.set_printoptions(edgeitems=3)
    np.core.arrayprint._line_width = 120


    print 'EFT coefficient matrix inversion'

    # Predefined values of Wilson coefs. for which the EFT tttt cross section was calculated.
    # One has to make sure that resulting matrix for sigma_i and sigma_ij is not degenerate
    c1 = [1, 0, 0, 0, 0]
    c2 = [0, 1, 0, 0, 0]
    c3 = [0, 0, 1, 0, 0]
    c4 = [0, 0, 0, 1, 0]
    c5 = [0, 0, 0, 0, 1]
    c6 = [1, 1, 1, 1, 1]
    c7 = [-1, -1, 1, 1, 1]
    c8 = [-1, -1, 1, 0, 1]
    c9 = [0, 1, 0, 0, -1]
    c10 = [0, 1, 1, 1, 0]
    c11 = [1, 0, -1, 1, 0]
    c12 = [-1, 0, 0, 1, -1]
    c13 = [-1, 0, 0, -1, 1]
    c14 = [0, 1, -1, 1, -1]
    c15 = [0, 1, 0, -1, 0]
    c16 = [0, 0, -1, -1, 1]
    c17 = [1, -1, 0, -1, 0]
    c18 = [1, 1, 0, -1, 1]
    c19 = [0, 1, 0, -1, 1]
    c20 = [1, -1, -1, 0, 1]

    # Fill matrix for sigma_i and sigma_ij
    A = np.array([gen_row(c1),gen_row(c2),gen_row(c3),gen_row(c4),gen_row(c5),gen_row(c6),gen_row(c7),gen_row(c8),gen_row(c9),gen_row(c10),gen_row(c11),gen_row(c12),gen_row(c13),gen_row(c14),gen_row(c15),gen_row(c16),gen_row(c17),gen_row(c18),gen_row(c19),gen_row(c20)])
    #print A

    # EFT cross section values for different values of Ci in the vectors c1, c2, c3, ... c19, c20
    MG_SM = [0.01557, 0.01564, 0.0102, 0.01116, 0.01022, 0.02873, 0.0203, 0.01704, 0.01527, 0.02066, 0.0168, 0.01741, 0.01567, 0.01342, 0.01839, 0.01208, 0.02113, 0.0283, 0.01983, 0.02386] #pb
    b = []

    # Subtract SM tttt cross section from EFT
    for i in range(0,len(MG_SM)): MG_SM[i] = MG_SM[i]*1000. - sig_SM
    b = np.array(MG_SM)
    #print b

    # solve linear system of equations for sigma_i, sigma_ij coefficients
    S = np.linalg.solve(A, b)
    #print S

    # solution sigma_i, sigma_ij
    sig_i = [S[0], S[1], S[2], S[3], S[4]];
    sig_ij = [S[5], S[10], S[11], S[12], S[13], S[10], S[6], S[14], S[15], S[16], S[11], S[14], S[7], S[17], S[18], S[12], S[15], S[17], S[8], S[19], S[13], S[16], S[18], S[19], S[9]];
    #print sig_i
    #print sig_ij

    #Plots with limit contours
    #make_plot1d(S)
    make_plot2d(S)
    
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
        yband = np.array([2.0*sig_SM, 2.0*sig_SM])
        band_half_w_1s_pos = np.array([(1.0)*sig_SM,(1.0)*sig_SM])
        band_half_w_1s_neg = np.array([(0.7)*sig_SM,(0.7)*sig_SM])
        band_half_w_2s_pos = np.array([(4.5)*sig_SM,(4.5)*sig_SM])
        band_half_w_2s_neg = np.array([(1.0)*sig_SM,(1.0)*sig_SM])

        plt.plot(xband, yband, 'k')
        plt.fill_between(xband, yband-band_half_w_2s_neg, yband+band_half_w_2s_pos, facecolor='green')
        plt.fill_between(xband, yband-band_half_w_1s_neg, yband+band_half_w_1s_pos, facecolor='yellow')

        plt.xlabel(r'$c_{R}$', labelpad=2)
        plt.ylabel(r'$\sigma_{t\bar{t}t\bar{t}}$ (fb)')

        two_s_band = mpatches.Patch(color='green', label=r'2 s.d.')
        one_s_band = mpatches.Patch(color='yellow', label=r'1 s.d.')
        limit_line = mlines.Line2D([], [], color='black', label='SL+OS+SS Combined')
        eft_line   = mlines.Line2D([], [], color='blue', label='EFT')
        plt.legend([eft_line,limit_line, one_s_band, two_s_band],['EFT','SL+OS+SS\nCombined','1 s.d.','2 s.d.'],loc=4,prop={'size':12})

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
        level_2015 = [5.9*sig_SM] 
        level_2016_sl = [2.0*sig_SM] 
        contour_2015 = plt.contour(X, Y, Z, level_2015, colors='r',linewidths=np.arange(3.9, 4, .5),linestyles='dashed')
        plt.clabel(contour_2015, colors = 'r', fmt = 'SL+OS', fontsize=12)
        contour_2016sl = plt.contour(X, Y, Z, level_2016_sl, colors='r',linewidths=np.arange(3.9, 4, .5))
        plt.clabel(contour_2016sl, colors = 'r', fmt = 'SL+OS+SS', fontsize=12)
        contour_filled = plt.contourf(X, Y, Z, 100)
        plt.colorbar(contour_filled)
        plt.title('')
        plt.xlabel(r'$C_{R}$')
        plt.ylabel(r'$C_{L}^{(1)}$')
        plt.show()
    
if __name__ == "__main__":
        main()
