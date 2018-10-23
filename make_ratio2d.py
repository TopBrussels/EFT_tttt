from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM
import numpy as np
import ROOT as rt

number_of_points_in_direction = [200,200]
max_range = [  4.0*np.pi,  4.0*np.pi ]
min_range = [ -4.0*np.pi, -4.0*np.pi ]
step_size = [ (max_range[0] - min_range[0])/number_of_points_in_direction[0], (max_range[1] - min_range[1])/number_of_points_in_direction[1] ]

rt.gStyle.SetOptStat(0)
hist2d = rt.TH2D('hist2d','fractional difference 13 TeV;C_{O_{R}};C_{O_{L1}}',number_of_points_in_direction[0],min_range[0]-step_size[0]/2.,max_range[0]+step_size[0]/2.,
                                                         number_of_points_in_direction[1],min_range[1]-step_size[1]/2.,max_range[1]+step_size[1]/2.)
eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
for i in range(0,number_of_points_in_direction[0]+1):
    for j in range(0,number_of_points_in_direction[0]+1):
        x = min_range[0]+i*step_size[0]
        y = min_range[1]+j*step_size[1]
        xs = eft.gen_eft_xs([x,y,0.,0.,0.])/sig_SM
        #print i,j,x,y
        hist2d.SetBinContent(i,j,xs)
#hist2d.Draw("colz text")
outfile = rt.TFile.Open("out.root","RECREATE")
outfile.cd()
hist2d.Write()
outfile.Close()
