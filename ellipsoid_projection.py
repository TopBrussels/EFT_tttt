import numpy as np

from mg_calculations import *
from eft_coefficients import EftPredictions

# see http://www.am.ub.edu/~robert/Documents/quadric.pdf

eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
S = eft.S

list_reference_wilson_coefs = [1.0]*5
np_reference_wilson_coefs = np.array(list_reference_wilson_coefs)

print "Reference cross section: ", eft.gen_eft_xs(list_reference_wilson_coefs) 

mat_result = np.dot( np_reference_wilson_coefs.T, np.dot(eft.Sigma2_matr, np_reference_wilson_coefs)) + 2.0*np.dot(eft.Sigma1_vec.T,np_reference_wilson_coefs) + eft.sig_sm
print "Matrix form result: ", mat_result

upper_limit_50 = 2.0312

Ak_matrices = []
Zk_vectors = []
Ck_vectors = []
Dk_vectors = []
Ek_vectors = []
p2,p1,p0 = [],[],[]

# independent_coef = 4
# eft.Sigma1_vec = eft.Sigma1_vec[np.ix_([independent_coef])]
# eft.Sigma2_matr = eft.Sigma2_matr[np.ix_([independent_coef],[independent_coef])]
# print eft.Sigma1_vec
# print eft.Sigma2_matr

N_wilsons = 5
for k in range(N_wilsons):
    Ak_mat = eft.Sigma2_matr.copy()
    Zk_vec = -eft.Sigma2_matr.copy()[:,k]
    Ck_vec = -eft.Sigma1_vec.copy()
    for i in range(N_wilsons):
        for j in range(N_wilsons):
            if i == k: Ak_mat[i,j] = 0
            if j == k: Ak_mat[i,j] = 0
    Ak_mat[k,k]=1
    Zk_vec[k]=1
    Ck_vec[k]=0
    Dk_vec = np.linalg.inv(Ak_mat).dot(Zk_vec)
    Ek_vec = np.linalg.inv(Ak_mat).dot(Ck_vec)
    # print "Ak:",Ak_mat
    # print "Zk:",Zk_vec
    # print "Ck:",Ck_vec
    # print "Dk:",Dk_vec
    # print "Ek:",Ek_vec
    Ak_matrices.append(Ak_mat)
    Zk_vectors.append(Zk_vec)
    Ck_vectors.append(Ck_vec)
    Dk_vectors.append(Dk_vec)
    Ek_vectors.append(Ek_vec)
    p2.append(eft.Sigma2_matr.copy()[k,:].dot(Dk_vec))
    p1.append(eft.Sigma2_matr.copy()[k,:].dot(Ek_vec) + eft.Sigma1_vec.copy().dot(Dk_vec) + eft.Sigma1_vec[k])
    p0.append(eft.Sigma1_vec.copy().dot(Ek_vec)+sig_SM)
# print p2
# print p1
# print p0

for k in range(N_wilsons):
    limits = np.roots([p2[k],p1[k],p0[k] - upper_limit_50*sig_SM])
    print "Marginal C{} range:  ".format(k), limits
