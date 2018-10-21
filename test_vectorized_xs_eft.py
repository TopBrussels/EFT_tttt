from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM
import numpy as np
import pprint as pp
number_of_points_in_direction = [20,20,20,20,20]
number_of_reduced_points_in_one_direction = 3
pointlist0 = np.linspace(-2.0, 2.0, num=number_of_points_in_direction[0])
pointlist1 = np.linspace(-9.0, 9.0, num=number_of_points_in_direction[1])
pointlist2 = np.linspace(-30.0, 30.0, num=number_of_points_in_direction[2])
pointlist3 = np.linspace(-3.0, 3.0, num=number_of_points_in_direction[3])
pointlist4 = np.linspace(-6.0, 6.0, num=number_of_points_in_direction[4])

# C0,C1,C2=np.meshgrid(pointlist0,pointlist1,pointlist2,sparse=False,indexing='ij')
# print C0
# print '-'*30
# print C0**2
# print '*'*30
# print C1
# print '*'*30
# print C2

upper_lim = 2.0312 * sig_SM
max_lim = {'C0':-999,'C1':-999,'C2':-999,'C3':-999,'C4':-999}
min_lim = {'C0':999,'C1':999,'C2':999,'C3':999,'C4':999}
C0,C1,C2,C3,C4 = np.meshgrid(pointlist0,pointlist1,pointlist2,pointlist3,pointlist4,sparse=False,indexing='ij')
eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
xs = eft.vgen_eft_xs(C0, C1, C2, C3, C4)
print "Number of cross section points: ", len(xs.ravel())
for i in range(number_of_points_in_direction[0]):
    for j in range(number_of_points_in_direction[1]):
        for k in range(number_of_points_in_direction[2]):
            for l in range(number_of_points_in_direction[3]):
                for m in range(number_of_points_in_direction[4]):
                    c0 = C0[i,j,k,l,m]
                    c1 = C1[i,j,k,l,m]
                    c2 = C2[i,j,k,l,m]
                    c3 = C3[i,j,k,l,m]
                    c4 = C4[i,j,k,l,m]
                    # ------------------
                    # c0 = 0.
                    # c1 = 0.
                    # c2 = 0.
                    # c3 = 0.
                    # c4 = 0.
                    if (xs[i,j,k,l,m] - upper_lim)<0: 
                        if max_lim['C0']<c0: max_lim['C0'] = c0
                        if min_lim['C0']>c0: min_lim['C0'] = c0
                        if max_lim['C1']<c1: max_lim['C1'] = c1
                        if min_lim['C1']>c1: min_lim['C1'] = c1
                        if max_lim['C2']<c2: max_lim['C2'] = c2
                        if min_lim['C2']>c2: min_lim['C2'] = c2
                        if max_lim['C3']<c3: max_lim['C3'] = c3
                        if min_lim['C3']>c3: min_lim['C3'] = c3
                        if max_lim['C4']<c4: max_lim['C4'] = c4
                        if min_lim['C4']>c4: min_lim['C4'] = c4
print "Marginal C0 range: ", min_lim['C0'], max_lim['C0']
print "Marginal C1 range: ", min_lim['C1'], max_lim['C1']
print "Marginal C2 range: ", min_lim['C2'], max_lim['C2']
print "Marginal C3 range: ", min_lim['C3'], max_lim['C3']
print "Marginal C4 range: ", min_lim['C4'], max_lim['C4']
# pp.pprint(C0)
# print '*'*30
# pp.pprint(C1)
# print '*'*30
# pp.pprint(C2)
# print '*'*30
# pp.pprint(C3)
# print '*'*30
# pp.pprint(C4)
# print '*'*30
# pp.pprint(C0*C1*C2*C3*C4)
