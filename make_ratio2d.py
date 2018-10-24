from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM, MG_SM_27TeV, sig_SM_27TeV, MG_SM_14TeV, sig_SM_14TeV
import numpy as np
import ROOT as rt
import CMS_lumi, tdrstyle
import itertools




tdrstyle.setTDRStyle()

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"


outfile = rt.TFile.Open("out.root","RECREATE")
outfile.cd()
wilsoncoeff_names=["C_{O_{R}}","C_{O_{L1}}","C_{O_{L8}}","C_{B1}","C_{B8}"]
thresholdval0=1.0
thresholdval1=0.30
thresholdval2=0.09
thresholdval3=0.01


number_of_points_in_direction = [100,100]
#max_range = [  4.0*np.pi,  4.0*np.pi ]
#min_range = [ -4.0*np.pi, -4.0*np.pi ]
scale=1.0

rt.gStyle.SetOptStat(0)


myfunc=rt.TF1("myfunc","[0]*x*x+[1]*x+[2]",-4.0*np.pi, 4.0*np.pi)

eft13 = EftPredictions(wilson_coefficients, MG_SM, sig_SM)
eft14 = EftPredictions(wilson_coefficients, MG_SM_14TeV, sig_SM_14TeV)
eft27 = EftPredictions(wilson_coefficients, MG_SM_27TeV, sig_SM_27TeV)

# 1D first

for iwilson in range(0,5):
    number_of_points_in_direction = 10000
    max_range =  4.0*np.pi
    min_range = -4.0*np.pi
    step_size =  (max_range - min_range)/number_of_points_in_direction
    histo1d13 = rt.TGraph(number_of_points_in_direction+1)
    histo1d13.SetName("13 TeV")
    histo1d14 = rt.TGraph(number_of_points_in_direction+1)
    histo1d14.SetName("14 TeV")
    histo1d27 = rt.TGraph(number_of_points_in_direction+1)
    histo1d27.SetName("27 TeV")
    for i in range(0,number_of_points_in_direction+1):
        x = min_range+i*step_size
        workarrayhorizontal=[0,0,0,0,0]
        workarrayhorizontal[iwilson]=x
        histo1d13.SetPoint(i,x,eft13.gen_eft_xs(workarrayhorizontal)/eft13.gen_eft_xs([0,0,0.,0.,0.]))
        histo1d14.SetPoint(i,x,eft14.gen_eft_xs(workarrayhorizontal)/eft14.gen_eft_xs([0,0,0.,0.,0.]))
        histo1d27.SetPoint(i,x,eft27.gen_eft_xs(workarrayhorizontal)/eft27.gen_eft_xs([0,0,0.,0.,0.]))

    canv = rt.TCanvas()
    histo1d13.Draw("la")
    histo1d14.Draw("l")
    histo1d27.Draw("l")
    canv.Update()
#    histo1d13.Fit(myfunc,"Q","goff")

    for histowork in histo1d13,histo1d14,histo1d27:
        histowork.Fit(myfunc,"Q","goff")
        for workval in thresholdval0,thresholdval1,thresholdval2,thresholdval3:
            a2=myfunc.GetParameter(0)
            b2=myfunc.GetParameter(1)
            c2=myfunc.GetParameter(2)-workval-1
            det2 =rt.TMath.Sqrt((b2*b2)-(4*a2*c2))
            print wilsoncoeff_names[iwilson], histowork.GetName(), workval, round((-b2+det2)/(2*a2),2), round((-b2-det2)/(2*a2),2)
#, round(histowork.Eval((-b2+det2)/(2*a2)),2),round(histowork.Eval((-b2-det2)/(2*a2)),2)

for combo in itertools.combinations(range(0,5),2) :
    print combo[0],combo[1],wilsoncoeff_names[combo[0]],wilsoncoeff_names[combo[1]]

    # do a simple loop over the diagonal to determine range first

    number_of_points_in_direction = [100,100]
    max_range = [  4.0*np.pi,  4.0*np.pi ]
    min_range = [ -4.0*np.pi, -4.0*np.pi ]

    
    step_size = [ (max_range[0] - min_range[0])/number_of_points_in_direction[0], (max_range[1] - min_range[1])/number_of_points_in_direction[1] ]
    
    histo1d1=rt.TGraph(number_of_points_in_direction[0]+1)
    histo1d2=rt.TGraph(number_of_points_in_direction[1]+1)
    # scan X and Y separately first:
    workrange=[0,0]
    for i in range(0,number_of_points_in_direction[0]+1):
        x = min_range[0]+i*step_size[0]
        workarrayhorizontal=[0,0,0,0,0]
        workarrayhorizontal[combo[0]]=x
        workarrayvertical=[0,0,0,0,0]
        workarrayvertical[combo[1]]=x
        
        histo1d1.SetPoint(i,x,eft13.gen_eft_xs(workarrayhorizontal)/eft13.gen_eft_xs([0,0,0.,0.,0.]))
        histo1d2.SetPoint(i,x,eft13.gen_eft_xs(workarrayvertical)/eft13.gen_eft_xs([0,0,0.,0.,0.]))
    

    canv = rt.TCanvas()
    histo1d1.Draw("la")
    histo1d2.Draw("l")
    histo1d1.GetHistogram().SetXTitle(wilsoncoeff_names[combo[0]]+","+wilsoncoeff_names[combo[1]])
    histo1d1.Fit(myfunc,"Q","goff")

    a1=myfunc.GetParameter(0)
    b1=myfunc.GetParameter(1)
    c1=myfunc.GetParameter(2)-thresholdval0-1
    det1=rt.TMath.Sqrt((b1*b1)-(4*a1*c1))

    histo1d2.Fit(myfunc,"Q","goff")
    a2=myfunc.GetParameter(0)
    b2=myfunc.GetParameter(1)
    c2=myfunc.GetParameter(2)-thresholdval0-1
    det2 =rt.TMath.Sqrt((b2*b2)-(4*a2*c2))
    # assuming a parabola, analytically solve range points:
    print wilsoncoeff_names[combo[0]], det1, (-b1+det1)/(2*a1), (-b1-det1)/(2*a1), histo1d1.Eval((-b1+det1)/(2*a1)),histo1d1.Eval((-b1-det1)/(2*a1))
    print wilsoncoeff_names[combo[1]], det2, (-b2+det2)/(2*a2), (-b2-det2)/(2*a2), histo1d2.Eval((-b2+det2)/(2*a2)),histo1d2.Eval((-b2-det2)/(2*a2))

    workrange[0]=max(abs((-b1+det1)/(2*a1)),abs((-b1-det1)/(2*a1)))
    workrange[1]=max(abs((-b2+det2)/(2*a2)),abs((-b2-det2)/(2*a2)))

    canv.Update()

    if workrange[0]<0.1:
        workrange[0]=0.1
    if workrange[1]<0.1:
        workrange[1]=0.1

    print "scanned range, setting to ",workrange
    max_range=[1.2*workrange[0],2*workrange[1]]
    min_range=[-1.2*workrange[0],-1.2*workrange[1]]
    number_of_points_in_direction[0]*=10
    number_of_points_in_direction[1]*=10

    step_size = [ (max_range[0] - min_range[0])/number_of_points_in_direction[0], (max_range[1] - min_range[1])/number_of_points_in_direction[1] ]
#    print max_range,min_range,step_size
    
    hist2d = rt.TH2D('hist2d_'+str(combo[0])+str(combo[1]),'fractional difference 13 TeV;'+wilsoncoeff_names[combo[0]]+';'+wilsoncoeff_names[combo[1]],number_of_points_in_direction[0],min_range[0]-step_size[0]/2.,max_range[0]+step_size[0]/2.,
                 number_of_points_in_direction[1],min_range[1]-step_size[1]/2.,max_range[1]+step_size[1]/2.)

    hist2d.SetTitle("")
    scenario0=hist2d.Clone("scenario0_"+str(combo[0])+str(combo[1]))
    scenario1=hist2d.Clone("scenario1_"+str(combo[0])+str(combo[1]))
    scenario2=hist2d.Clone("scenario2_"+str(combo[0])+str(combo[1]))
    scenario3=hist2d.Clone("scenario3_"+str(combo[0])+str(combo[1]))



    for i in range(0,number_of_points_in_direction[0]+1):
        for j in range(0,number_of_points_in_direction[0]+1):
            x = min_range[0]+i*step_size[0]
            y = min_range[1]+j*step_size[1]
            workarray=[0,0,0,0,0]
            workarray[combo[0]]=x
            workarray[combo[1]]=y
            #            print workarray
            #        print i,j,x,y,xs
            hist2d.SetBinContent(i,j,eft13.gen_eft_xs(workarray)/eft13.gen_eft_xs([0,0,0.,0.,0.]))
            if eft13.gen_eft_xs(workarray)/eft13.gen_eft_xs([0,0,0.,0.,0.]) < 1.+thresholdval0 :
                #            print i,j,x,y,xs
                scenario0.SetBinContent(i,j,0.)
            else:
                scenario0.SetBinContent(i,j,1.)
            if eft14.gen_eft_xs(workarray)/eft14.gen_eft_xs([0,0,0.,0.,0.]) < 1.+thresholdval1 :
                #            print i,j,x,y,xs
                scenario1.SetBinContent(i,j,0.)
            else:
                scenario1.SetBinContent(i,j,1.)
            if eft14.gen_eft_xs(workarray)/eft14.gen_eft_xs([0,0,0.,0.,0.]) < 1.+thresholdval2 :
                #            print i,j,x,y,xs
                scenario2.SetBinContent(i,j,0.)
            else:
                scenario2.SetBinContent(i,j,1.)
            if eft27.gen_eft_xs(workarray)/eft27.gen_eft_xs([0,0,0.,0.,0.]) < 1.+thresholdval3 :
                #            print i,j,x,y,xs
                scenario3.SetBinContent(i,j,0.)
            else:
                scenario3.SetBinContent(i,j,1.)

    legendstring0 = "LHC 13TeV 35.9 fb^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval0,0))+"%"
    legendstring1 = "HL-LHC 3 ab^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval1,0))+"%"
    legendstring2 = "HL-LHC 0.3 ab^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval2,0))+"%"
    legendstring3 = "HE-LHC 15 ab^{-1}, #Delta#sigma_{tttt} : "+str(round(100*thresholdval3,0))+"%"

    hist2d.SetTitle("")
    #hist2d.Draw("CONT1Z")
    #rt.gPad.SetLogz()
    hist2d.Draw("axis")
    scenario0.SetLineColor(rt.kGrey)
    scenario0.SetLineStyle(3)
    scenario0.Draw("CONT3SAME")
    scenario1.SetLineColor(rt.kAzure+1)
    scenario1.Draw("CONT3SAME")
    scenario2.SetLineColor(rt.kGreen+3)
    scenario2.Draw("CONT3SAME")
    scenario3.SetLineColor(rt.kRed)
    scenario3.Draw("CONT3SAME")
    rt.gPad.Update()
    leg = rt.TLegend(0.5,0.75,0.93,0.93)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(scenario0,legendstring0,"l")
    leg.AddEntry(scenario1,legendstring1,"l")
    leg.AddEntry(scenario2,legendstring2,"l")
    leg.AddEntry(scenario3,legendstring3,"l")
    leg.Draw("same")
    CMS_lumi.CMS_lumi(rt.gPad,0,11)
    rt.gPad.Update()
    rt.gPad.Print("plot"+str(combo[0])+str(combo[1])+".png")
    rt.gPad.Print("plot"+str(combo[0])+str(combo[1])+".pdf")
    rt.gPad.Print("plot"+str(combo[0])+str(combo[1])+".C")
    outfile.cd()
hist2d.Write()
outfile.Close()

