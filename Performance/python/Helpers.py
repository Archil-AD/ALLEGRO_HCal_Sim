#!/usr/bin/python
import os
import math
import time
import ROOT as root
from ROOT import TFile, TTree, TPad, TH1
from array import array
import numpy as np

def plot(hist, opt):
    if not os.path.exists('plots'):
      os.makedirs('plots')
    ######################
    Logy = False
    Logx = False
    ymin = 0.0001
    ymax = 1
    labels = []
    name = ''
    ######################
    if hasattr(opt,'Name'):
      name = opt.Name

    if hasattr(opt,'Labels'):
      labels = opt.Labels

    if hasattr(opt,'Logy'):
      Logy = opt.Logy
    if hasattr(opt,'Logx'):
      Logx = opt.Logx

    if hasattr(opt,'xTitle'):
      hist.GetXaxis().SetTitle(opt.xTitle);

    if hasattr(opt,'yTitle'):
      hist.GetYaxis().SetTitle(opt.yTitle);

    if hasattr(opt, 'Normalize'):
      if opt.Normalize:
        hist.Scale(1./hist.Integral(0,-1))
        hist.GetYaxis().SetTitle('Normalized entries')

    if hasattr(opt,'yRange'):
      ymin = opt.yRange[0]
      ymax = opt.yRange[1]
    else:
      ymax = hist.GetMaximum()*1.2

    if not Logx and hasattr(opt,'xRange'):
      hist.GetXaxis().SetRangeUser(opt.xRange[0],opt.xRange[1]);

    if not Logy:
      hist.GetYaxis().SetRangeUser(ymin,ymax);
    else:
      hist.SetMinimum(ymin)
      hist.SetMaximum(ymax)
    ###########################

    c1 = root.TCanvas("c1","",0,0,800,600)
    root.gStyle.SetLegendBorderSize(0)

    if Logx:
      root.gPad.SetLogx()
    if Logy:
      root.gPad.SetLogy()

    hist.GetXaxis().SetLabelSize(0.055);
    hist.GetXaxis().SetTickLength(0.02);
    hist.GetXaxis().SetTitleSize(0.055);
    if Logx:
      hist.GetXaxis().SetMoreLogLabels()

    hist.GetYaxis().SetTickLength(0.02);
    hist.GetYaxis().SetTitleOffset(1.2);
    hist.GetYaxis().SetTitleSize(0.055);
    hist.GetYaxis().SetLabelOffset(0.005);
    hist.GetYaxis().SetLabelSize(0.055);

    if hasattr(opt, 'DrawOpt'):
      hist.Draw(opt.DrawOpt)
    else: hist.Draw("E1X0")

    if hasattr(opt, 'Hists'):
      hists = opt.Hists
      markers = [20,21,22,23]
      colors = [1,Color(0),Color(1),Color(2)]
      leg = root.TLegend(0.18,0.75,0.38,0.95);
      leg.SetHeader("");
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetShadowColor(0);
      leg.SetLineColor(0);
      leg.SetLineWidth(0);
      leg.SetTextSize(0.045);
      leg.SetTextFont(42);

      for i in range(0,len(hists)):
        hists[i].SetLineColor(colors[i])
        hists[i].SetMarkerColor(colors[i])
        hists[i].SetMarkerStyle(markers[i])
        hists[i].SetMarkerSize(1.5)
        hists[i].Draw("same hist")
        leg.AddEntry(hists[i],hists[i].GetTitle(),"f");
      leg.Draw();

    ALLEGRO_LABEL(0.56,0.85,"Simulation",0.2)
    for i in range(len(labels)):
      myText( 0.56,0.842-(i+1)*0.055,1,"#scale[1]{%s}" % labels[i]);

    if hasattr(opt, 'Function'):
      if opt.Function:
        opt.Function.SetLineColor(Color(0))
        opt.Function.Draw("same")

    c1.Print(f'plots/{name}')


def plotGraph(graphs, opt):
    if not os.path.exists('plots'):
      os.makedirs('plots')
    hist = graphs[0]
    ######################
    Logy = False
    Logx = False
    ymin = 0.0001
    ymax = 1
    xmin = 0
    xmax = 150
    labels = []
    name = ''
    ######################
    if hasattr(opt,'Name'):
      name = opt.Name

    if hasattr(opt,'Labels'):
      labels = opt.Labels

    if hasattr(opt,'Logy'):
      Logy = opt.Logy
    if hasattr(opt,'Logx'):
      Logx = opt.Logx

    if hasattr(opt,'xTitle'):
      hist.GetXaxis().SetTitle(opt.xTitle);

    if hasattr(opt,'yTitle'):
      hist.GetYaxis().SetTitle(opt.yTitle);

    if hasattr(opt,'yRange'):
      ymin = opt.yRange[0]
      ymax = opt.yRange[1]
    else:
      ymax = hist.GetMaximum()*1.2

    if hasattr(opt,'xRange'):
      xmin = opt.xRange[0]
      xmax = opt.xRange[1]
      hist.GetXaxis().SetLimits(xmin,xmax)

    hist.SetMinimum(ymin)
    hist.SetMaximum(ymax)
    ###########################

    c1 = root.TCanvas("c1","",0,0,800,600)
    root.gStyle.SetLegendBorderSize(0)

    if Logx:
      root.gPad.SetLogx()
    if Logy:
      root.gPad.SetLogy()

    hist.GetXaxis().SetLabelSize(0.055);
    hist.GetXaxis().SetTickLength(0.02);
    hist.GetXaxis().SetTitleSize(0.055);
    if Logx:
      hist.GetXaxis().SetMoreLogLabels()

    hist.GetYaxis().SetTickLength(0.02);
    hist.GetYaxis().SetTitleOffset(1.2);
    hist.GetYaxis().SetTitleSize(0.055);
    hist.GetYaxis().SetLabelOffset(0.005);
    hist.GetYaxis().SetLabelSize(0.055);
    hist.SetMarkerSize(1.5);

    hist.Draw("AP")
    ALLEGRO_LABEL(0.56,0.85,"Simulation",0.2)
    for i in range(len(labels)):
      myText( 0.56,0.842-(i+1)*0.055,1,"#scale[1]{%s}" % labels[i]);

    if len(graphs) > 1 and len(graphs) < 5:
      markers = [20,21,22,23]
      colors = [1,Color(0),Color(1),Color(2)]
      leg = root.TLegend(0.18,0.75,0.38,0.95);
      leg.SetHeader("");
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetShadowColor(0);
      leg.SetLineColor(0);
      leg.SetLineWidth(0);
      leg.SetTextSize(0.045);
      leg.SetTextFont(42);
      leg.AddEntry(graphs[0],graphs[0].GetTitle(),"p");

      for i in range(1,len(graphs)):
        graphs[i].SetLineColor(colors[i])
        graphs[i].SetMarkerColor(colors[i])
        graphs[i].SetMarkerStyle(markers[i])
        graphs[i].SetMarkerSize(1.5)
        graphs[i].Draw("same p")
        leg.AddEntry(graphs[i],graphs[i].GetTitle(),"p");
      leg.Draw();

    line = root.TLine(xmin,1.0,xmax,1.)
    line.SetLineColor(root.kBlack);
    line.SetLineStyle(7);
    line.SetLineWidth(2);
    line.Draw("same");

    if hasattr(opt, 'Function'):
      opt.Function[0].SetLineColor(Color(0))
      opt.Function[0].Draw("same")
      leg = root.TLegend(0.18,0.8,0.38,0.95);
      leg.SetHeader("");
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetShadowColor(0);
      leg.SetLineColor(0);
      leg.SetLineWidth(0);
      leg.SetTextSize(0.045);
      leg.SetTextFont(42);
      leg.AddEntry(opt.Function[0],opt.Function[1],"l");
      leg.Draw()

    c1.Print(f'plots/{name}')



def plotCompare(hist0 ,hist1, opt, hist2 = [None,None]):
    if not os.path.exists('plots'):
      os.makedirs('plots')
    h0 = hist0[0]
    h1 = hist1[0]
    h2 = hist2[0]
    ######################
    Logy = False
    Logx = False
    rymin = 0.85
    rymax = 1.15
    ymin = 0.0001
    ymax = 1
    labels = []
    name = ''
    ######################
    if hasattr(opt,'Name'):
      name = opt.Name

    if hasattr(opt,'Labels'):
      labels = opt.Labels

    if hasattr(opt,'Logy'):
      Logy = opt.Logy
    if hasattr(opt,'Logx'):
      Logx = opt.Logx

    if hasattr(opt,'RatioRange'):
      rymin = opt.RatioRange[0]
      rymax = opt.RatioRange[1]

    if hasattr(opt,'xTitle'):
      h0.GetXaxis().SetTitle(opt.xTitle);

    if hasattr(opt,'yTitle'):
      h0.GetYaxis().SetTitle(opt.yTitle);

    if hasattr(opt, 'Normalize'):
      if opt.Normalize:
        h0.Scale(1./h0.Integral(0,-1))
        h1.Scale(1./h1.Integral(0,-1))
        if h2:
          h2.Scale(1./h2.Integral(0,-1))
        h0.GetYaxis().SetTitle('Normalized entries')

    if hasattr(opt,'Ndivisions'):
      h0.GetXaxis().SetNdivisions(opt.Ndivisions)

    if hasattr(opt,'Rebin'):
      h0.Rebin(opt.Rebin)
      h1.Rebin(opt.Rebin)

    if hasattr(opt,'yRange'):
      ymin = opt.yRange[0]
      ymax = opt.yRange[1]
    else:
      max0 = h0.GetMaximum()
      max1 = h1.GetMaximum()
      if max0 > max1:
        ymax= max0*1.5
      else:
        ymax= max1*1.5

    if not Logx and hasattr(opt,'xRange'):
      h0.GetXaxis().SetRangeUser(opt.xRange[0],opt.xRange[1]);
    ###########################


    root.gStyle.SetLegendBorderSize(0)
    c1 = root.TCanvas("c1","multipads",700,700)
    pad1 = root.TPad("pad1","pad1",0.,0.3,1.,1.);
    pad2 = root.TPad("pad2","pad2",0.,0.,1.,0.35);

    pad1.SetTopMargin(0.075);
    pad2.SetTopMargin(0.00);
    pad1.SetBottomMargin(0.085);
    pad2.SetBottomMargin(0.34);

    pad2.SetGridy();
    pad2.SetFillColor(0);
    #pad2.SetGridy();
    pad2.SetTickx(1);
    pad2.SetTicky(1);
    if Logx:
      pad1.SetLogx()
      pad2.SetLogx()
    if Logy: 
      pad1.SetLogy()

    pad1.SetLeftMargin(0.15);
    pad2.SetLeftMargin(0.15);

    pad1.SetRightMargin(0.075);
    pad2.SetRightMargin(0.075);

    pad1.SetFrameBorderMode(0);
    pad2.SetFrameBorderMode(0);

    pad1.SetFillStyle(0);
    pad2.SetFillStyle(0);

    pad1.Draw();
    pad2.Draw();
    pad1.cd();

    h0.SetLineStyle(1);
    h0.SetLineColor(1);
    h0.SetMarkerStyle(20)
    h0.SetMarkerSize(1.5);
    h0.SetMarkerColor(1);
    h1.SetLineStyle(Color(0));
    h1.SetLineColor(Color(0));
    h1.SetMarkerStyle(22)
    h1.SetMarkerSize(1.5);
    h1.SetMarkerColor(Color(0));
    if h2:
      h2.SetLineStyle(1);
      h2.SetLineColor(Color(6));
      h2.SetMarkerStyle(23)
      h2.SetMarkerColor(Color(6));
      h2.SetMarkerSize(1.5);

    h0.GetXaxis().SetLabelOffset(0.0);
    h0.GetXaxis().SetLabelSize(0.0);
    h0.GetXaxis().SetTickLength(0.025);
    h0.GetXaxis().SetTitleSize(0.);

    h0.GetYaxis().SetTickLength(0.015);
    h0.GetYaxis().SetTitleOffset(1);
    h0.GetYaxis().SetTitleSize(0.06);
    h0.GetYaxis().SetLabelOffset(0.005);
    h0.GetYaxis().SetLabelSize(0.055);

    if Logx:
      h0.GetXaxis().SetMoreLogLabels()

    if not Logy:
      h0.GetYaxis().SetRangeUser(ymin,ymax);
    else:
      h0.SetMinimum(ymin)
      h0.SetMaximum(ymax)

    h0.Draw();
    h1.Draw("same");
    if h2:
      h2.Draw("same");


    ALLEGRO_LABEL(0.59,0.85,"Internal",0.2)
    for i in range(len(labels)):
      myText( 0.59,0.842-(i+1)*0.055,1,"#scale[1]{%s}" % labels[i]);

    leg = root.TLegend(0.18,0.67,0.38,0.92);
    leg.SetHeader("");
    leg.AddEntry(h0, hist0[1],"p");
    leg.AddEntry(0, '',"");
    leg.AddEntry(h1, hist1[1],"p");
    if h2:
      leg.AddEntry(h2, hist2[1],"p");
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetShadowColor(0);
    leg.SetLineColor(0);
    leg.SetLineWidth(0);
    leg.SetTextSize(0.05);
    leg.SetTextFont(42);
    leg.Draw();

    ############
    pad2.cd()

    r0 = h0.Clone("r0")
    r0.SetLineColor(1)
    r0.Divide(h1)
#    for b in range(1,r0.GetNbinsX()+1):
#      if r0.GetBinContent(b) > 0:
#        r0.SetBinContent(b,1./r0.GetBinContent(b))
    if h2:
      r0.SetMarkerStyle(h1.GetMarkerStyle())
      r0.SetLineColor(h1.GetLineColor())
      r0.SetMarkerColor(h1.GetMarkerColor())
      r1 = h0.Clone("r1")
      r1.Divide(h2)
      r1.SetMarkerStyle(h2.GetMarkerStyle())
      r1.SetLineColor(h2.GetLineColor())
      r1.SetMarkerColor(h2.GetMarkerColor())
      for b in range(1,r1.GetNbinsX()+1):
        if r1.GetBinContent(b) > 0:
          r1.SetBinContent(b,1./r1.GetBinContent(b))


    r0.GetXaxis().SetTitle(h0.GetXaxis().GetTitle());
    r0.GetXaxis().SetTickLength(0.05);
    r0.GetYaxis().SetTickLength(0.02);

    r0.GetXaxis().SetTitleOffset(1.3);
    r0.GetXaxis().SetTitleSize(0.12);
    r0.GetXaxis().SetLabelOffset(0.025);
    r0.GetXaxis().SetLabelSize(0.11);
    r0.GetXaxis().LabelsOption("u");

    r0.GetYaxis().SetTitleSize(0.11);
    r0.GetYaxis().SetTitleOffset(0.55);
    r0.GetYaxis().SetLabelSize(0.11);
    r0.GetYaxis().SetLabelOffset(0.005);

#    r0.GetYaxis().SetTitle("Data22 / Data23");
    if not h2:
      r0.GetYaxis().SetTitle("Demo. / Legacy");
    else:
      r0.GetYaxis().SetTitle("Ratio to Data22");
    r0.GetYaxis().SetRangeUser(rymin,rymax);
    r0.GetYaxis().SetNdivisions(505);

    if Logx:
      r0.GetXaxis().SetLabelOffset(-0.01);

    xlow = h0.GetXaxis().GetBinLowEdge(1)
    xhigh = h0.GetXaxis().GetBinUpEdge(h0.GetNbinsX())
    if not Logx and hasattr(opt,'xRange'):
      xlow = opt.xRange[0]
      xhigh = opt.xRange[1]
      r0.GetXaxis().SetRangeUser(opt.xRange[0],opt.xRange[1]);

    r0.Draw()
    if h2:
      r1.Draw("same")

    if xlow < r0.GetBinLowEdge(1):
      xlow = r0.GetBinLowEdge(1)
    line = root.TLine(xlow,1.0,xhigh,1.)
    line.SetLineColor(root.kBlack);
    line.SetLineStyle(7);
    line.SetLineWidth(2);
    line.Draw("same");

    pad1.RedrawAxis()
    pad2.RedrawAxis()
    c1.RedrawAxis();
    c1.Print('plots/'+name)
    c1.Print('plots/'+name[:-4]+'.eps')
    c1.Print('plots/'+name[:-4]+'.C')
    c1.Print('plots/'+name[:-4]+'.root')


def IQnR(h):
    xq = array('d')  # position where to compute the quantiles in [0,1]
    xq.append(0.16)
    xq.append(0.5)
    xq.append(0.84)
    yq = array('d') # array to contain the quantiles
    yq.append(0.)
    yq.append(0.)
    yq.append(0.)

    h = OptimalRebin(h)

    h.GetQuantiles(3,yq,xq)

    return (yq[2]-yq[0])/(2.*yq[1])

def Median(h):
    xq = array('d')  # position where to compute the quantiles in [0,1]
    xq.append(0.16)
    xq.append(0.5)
    xq.append(0.84)
    yq = array('d') # array to contain the quantiles
    yq.append(0.)
    yq.append(0.)
    yq.append(0.)

    h = OptimalRebin(h)

    h.GetQuantiles(3,yq,xq)

    return yq[1]


def FitGaus(h, nsigma, mean, sigma):
    if h.GetEffectiveEntries() < 100 :
      return

    LIMIT=0.01;
    NITER=2;

    mean[0] =  h.GetMean()
    sigma[0] = h.GetRMS()

    h.Fit("gaus","QIME","", mean[0]-nsigma[0]*sigma[0], mean[0]+nsigma[1]*sigma[0])

    mean[0] =  h.GetFunction("gaus").GetParameter(1);
    sigma[0] = h.GetFunction("gaus").GetParameter(2);
    mean[1] =  h.GetFunction("gaus").GetParError(1);
    sigma[1] = h.GetFunction("gaus").GetParError(2);

    for i in range(NITER) :
      h.Fit("gaus","QIME","", mean[0]-nsigma[0]*sigma[0], mean[0]+nsigma[1]*sigma[0])
      if abs((h.GetFunction("gaus").GetParameter(1)-mean[0])/mean[0]) < LIMIT and abs((h.GetFunction("gaus").GetParameter(2)-sigma[0])/sigma[0]) < LIMIT :
        break
      else:
        mean[0] =  h.GetFunction("gaus").GetParameter(1)
        sigma[0] = h.GetFunction("gaus").GetParameter(2)
        mean[1] =  h.GetFunction("gaus").GetParError(1)
        sigma[1] = h.GetFunction("gaus").GetParError(2)

    return


def OptimalRebin( h, verbose = False ):
    method = 1
    N = h.GetEffectiveEntries();
    optWidth = 3.5*h.GetRMS()/root.TMath.Power(N,1.0/3);
    Nbins=int(h.GetNbinsX());
    fitRange = h.GetBinLowEdge(Nbins+1) - h.GetBinLowEdge(1);
    rebin=1;
    prevWidth=fitRange/Nbins;
    for i in range(1, Nbins) :
      if (Nbins%i!=0):
        continue;
      binWidth=fitRange/Nbins*i;

      if (method==1):
        if (binWidth<optWidth):
           rebin=i;
      elif (method==2):
        if (ROOT.TMath.Abs(binWidth-optWidth) < ROOT.TMath.Abs(prevWidth-optWidth)):
          rebin=i;
      else:
        rebin=i; # method 3

      if (binWidth>optWidth):
         break;
      prevWidth=binWidth;
    h.Rebin(rebin)
    if verbose:
      print ("\n%s\n  RMS: %.3f, Neff: %.3f\n" % (h.GetName(),h.GetRMS(),N) )
      print ("  Opt width: %6.3f, histo binwidth: %6.3f => Rebin: %d\n" % (optWidth,fitRange/Nbins,rebin) )
    return h

def Color(sequence_number):
  color_sequence = ["#3f90da","#ffa90e","#bd1f01","#94a4a2","#832db6","#a96b59","#e76300","#b9ac70","#717581","#92dadd"];
  return root.TColor.GetColor( color_sequence[ sequence_number%10 ] )

def ALLEGRO_LABEL(x, y, label = "Simulation", offset = 0.08):
    l = root.TLatex()  #l.SetTextAlign(12); l.SetTextSize(tsize);
    l.SetNDC()
    l.SetTextFont(62)
    #l.SetTextSize(0.042*size)
    l.DrawLatex(x,y,"ALLEGRO")
    l2 = root.TLatex()  #l.SetTextAlign(12);
    #l2.SetTextSize(0.042*size);
    l2.SetNDC()
    l2.DrawLatex(x+offset,y,label)

def myText(x, y, size, text):
    l = root.TLatex()  #l.SetTextAlign(12); l.SetTextSize(tsize);
    l.SetNDC()
    l.SetTextSize(0.04*size)
    l.SetTextColor(1)
    l.DrawLatex(x,y,text)
