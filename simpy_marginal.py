import sympy as sp
from mg_calculations import *
from eft_coefficients import EftPredictions

eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
S = eft.S
print(S)

sig_sm_tttt,C0,C1,C2,C3,C4 = sp.symbols("sigma_sm_tttt C0 C1 C2 C3 C4")
S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19 = sp.symbols("S0 S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19")

xs_eq = sig_sm_tttt +  C0 * S0 +  C1 * S1 +  C2 * S2 +  C3 * S3 +  C4 * S4 + \
               (C0 ** 2.) * S5 +  (C1 ** 2.) * S6 +  (C2 ** 2.) * S7 +  (C3 ** 2.) * S8 +  (C4 ** 2.) * S9 + \
               2. * C0 * C1 * S10 +  2. * C0 * C2 * S11 +  2. * C0 * C3 * S12 +  2. * C0 * C4 * S13 + \
               2. * C1 * C2 * S14 +  2. * C1 * C3 * S15 +  2. * C1 * C4 * S16 + \
               2. * C2 * C3 * S17 +  2. * C2 * C4 * S18 + \
               2. * C3 * C4 * S19
# print(xs_eq)
print('Analytic equation after substitution:',xs_eq.evalf(subs={sig_sm_tttt:sig_SM,C0:2.,C1:2.,C2:2.,C3:2.,C4:2.,
S0:S[0],
S1:S[1],
S2:S[2],
S3:S[3],
S4:S[4],
S5:S[5],
S6:S[6],
S7:S[7],
S8:S[8],
S9:S[9],
S10:S[10],
S11:S[11],
S12:S[12],
S13:S[13],
S14:S[14],
S15:S[15],
S16:S[16],
S17:S[17],
S18:S[18],
S19:S[19]
}))
print('Numerical correct equation:          ',eft.gen_eft_xs([2.,2.,2.,2.,2.]))
C0_ana_sol = sp.solve(xs_eq,C0)
print(C0_ana_sol)

# #First solution
d1C1 = sp.diff(C0_ana_sol[0],C1,1)
d1C1_subs = d1C1.subs({sig_sm_tttt:sig_SM,
S0:S[0],
S1:S[1],
S2:S[2],
S3:S[3],
S4:S[4],
S5:S[5],
S6:S[6],
S7:S[7],
S8:S[8],
S9:S[9],
S10:S[10],
S11:S[11],
S12:S[12],
S13:S[13],
S14:S[14],
S15:S[15],
S16:S[16],
S17:S[17],
S18:S[18],
S19:S[19]
})
d1C2 = sp.diff(C0_ana_sol[0],C2,1)
d1C2_subs = d1C2.subs({sig_sm_tttt:sig_SM,
S0:S[0],
S1:S[1],
S2:S[2],
S3:S[3],
S4:S[4],
S5:S[5],
S6:S[6],
S7:S[7],
S8:S[8],
S9:S[9],
S10:S[10],
S11:S[11],
S12:S[12],
S13:S[13],
S14:S[14],
S15:S[15],
S16:S[16],
S17:S[17],
S18:S[18],
S19:S[19]
})
d1C3 = sp.diff(C0_ana_sol[0],C3,1)
d1C3_subs = d1C3.subs({sig_sm_tttt:sig_SM,
S0:S[0],
S1:S[1],
S2:S[2],
S3:S[3],
S4:S[4],
S5:S[5],
S6:S[6],
S7:S[7],
S8:S[8],
S9:S[9],
S10:S[10],
S11:S[11],
S12:S[12],
S13:S[13],
S14:S[14],
S15:S[15],
S16:S[16],
S17:S[17],
S18:S[18],
S19:S[19]
})
d1C4 = sp.diff(C0_ana_sol[0],C4,1)
d1C4_subs = d1C4.subs({sig_sm_tttt:sig_SM,
S0:S[0],
S1:S[1],
S2:S[2],
S3:S[3],
S4:S[4],
S5:S[5],
S6:S[6],
S7:S[7],
S8:S[8],
S9:S[9],
S10:S[10],
S11:S[11],
S12:S[12],
S13:S[13],
S14:S[14],
S15:S[15],
S16:S[16],
S17:S[17],
S18:S[18],
S19:S[19]
})
print("First sol.:", d1C1,'\n\n', d1C2,'\n\n',d1C3,'\n\n',d1C4,'\n\n')
print("First sol. (after substitutions):", d1C1_subs,'\n\n', d1C2_subs,'\n\n',d1C3_subs,'\n\n',d1C4_subs,'\n\n')
print("Start numerical solving")
C0_1_sol = sp.nsolve((d1C1_subs, d1C2_subs, d1C3_subs, d1C4_subs), (C1, C2, C3, C4), (5.6, 12.5, 3., 12.5),prec=3,tol=1.0e-3)
print(C0_1_sol)
print(C0_ana_sol[0].evalf(subs={C1:C0_1_sol[0],C2:C0_1_sol[1],C3:C0_1_sol[2],C4:C0_1_sol[3],
sig_sm_tttt:sig_SM,
S0:S[0],
S1:S[1],
S2:S[2],
S3:S[3],
S4:S[4],
S5:S[5],
S6:S[6],
S7:S[7],
S8:S[8],
S9:S[9],
S10:S[10],
S11:S[11],
S12:S[12],
S13:S[13],
S14:S[14],
S15:S[15],
S16:S[16],
S17:S[17],
S18:S[18],
S19:S[19]
}))
#Second solution
# d2C1 = sp.diff(C0_ana_sol[1],C1,1)
# d2C2 = sp.diff(C0_ana_sol[1],C2,1)
# d2C3 = sp.diff(C0_ana_sol[1],C3,1)
# d2C4 = sp.diff(C0_ana_sol[1],C4,1)
# print("Second sol.:", d2C1,'\n\n',d2C2,'\n\n',d2C3,'\n\n',d2C4,'\n\n')
