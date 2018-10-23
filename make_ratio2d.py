from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM
import numpy as np
import ROOT as rt
import CMS_lumi, tdrstyle

tdrstyle.setTDRStyle()

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

thresholdval1=0.10
thresholdval2=0.03
thresholdval3=0.32
number_of_points_in_direction = [200,200]
#max_range = [  4.0*np.pi,  4.0*np.pi ]
#min_range = [ -4.0*np.pi, -4.0*np.pi ]
scale=0.1
max_range = [  scale*4.0*np.pi,  scale*4.0*np.pi ]
min_range = [ scale*-4.0*np.pi, scale*-4.0*np.pi ]


step_size = [ (max_range[0] - min_range[0])/number_of_points_in_direction[0], (max_range[1] - min_range[1])/number_of_points_in_direction[1] ]

rt.gStyle.SetOptStat(0)
hist2d = rt.TH2D('hist2d','fractional difference 13 TeV;C_{O_{R}};C_{O_{L1}}',number_of_points_in_direction[0],min_range[0]-step_size[0]/2.,max_range[0]+step_size[0]/2.,
                 number_of_points_in_direction[1],min_range[1]-step_size[1]/2.,max_range[1]+step_size[1]/2.)

hist2d.SetTitle("")
scenario1=hist2d.Clone("scenario1")
scenario2=hist2d.Clone("scenario2")
scenario3=hist2d.Clone("scenario3")

eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
for i in range(0,number_of_points_in_direction[0]+1):
    for j in range(0,number_of_points_in_direction[0]+1):
        x = min_range[0]+i*step_size[0]
        y = min_range[1]+j*step_size[1]
        xs = eft.gen_eft_xs([x,y,0.,0.,0.])/sig_SM
        #        print i,j,x,y,xs
        hist2d.SetBinContent(i,j,xs)
        if xs < 1.+thresholdval1 :
            #            print i,j,x,y,xs
            scenario1.SetBinContent(i,j,0.)
        else:
            scenario1.SetBinContent(i,j,1.)
        if xs < 1.+thresholdval2 :
            #            print i,j,x,y,xs
            scenario2.SetBinContent(i,j,0.)
        else:
            scenario2.SetBinContent(i,j,1.)
        if xs < 1.+thresholdval3 :
            #            print i,j,x,y,xs
            scenario3.SetBinContent(i,j,0.)
        else:
            scenario3.SetBinContent(i,j,1.)

legendstring1 = "HL-LHC 3 ab^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval1,0))+"%"
legendstring2 = "HE-LHC 15 ab^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval2,0))+"%"
legendstring3 = "HL-LHC 0.3 ab^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval3,0))+"%"

hist2d.SetTitle("")
#hist2d.Draw("CONT1Z")
#rt.gPad.SetLogz()
#hist2d.Draw("CONT1Z")
scenario1.SetLineColor(rt.kBlack)
scenario1.Draw("CONT3")
scenario2.SetLineColor(rt.kAzure+1)
scenario2.Draw("CONT3SAME")
scenario3.SetLineColor(rt.kRed)
scenario3.Draw("CONT3SAME")
rt.gPad.Update()
leg = rt.TLegend(0.5,0.75,0.95,0.95)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(scenario1,legendstring1,"l")
leg.AddEntry(scenario2,legendstring2,"l")
leg.AddEntry(scenario3,legendstring3,"l")
leg.Draw("same")
CMS_lumi.CMS_lumi(rt.gPad,0,11)
rt.gPad.Update()
rt.gPad.Print("plot.png")
outfile = rt.TFile.Open("out.root","RECREATE")
outfile.cd()
hist2d.Write()
outfile.Close()


