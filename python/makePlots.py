#!/usr/bin/python
import argparse
import math
import time
import sys
import ROOT as root
from ROOT import TFile, TTree, TH1F, TProfile
from array import array
import uncertainties as unc
import numpy as np
from scipy import stats
import os
import Helpers
import ctypes
class option:
    def __init__(self):
      self.name = ''

if __name__ == "__main__":
  root.gROOT.SetBatch(True) 
  # Activate the storage of the sum of squares of errors by default.
  root.TH1.SetDefaultSumw2(True)
  # Activate use of underflows and overflows in `Fill()` in the
  # computation of statistics (mean value, RMS) by default.
  root.TH1.StatOverflows(True)

  root.gROOT.SetStyle('ATLAS')


  parser = argparse.ArgumentParser(description='Process some integers.')
  parser.add_argument("--inputDir", dest='input', default="root://128.141.173.81:1094//home/data/FCC/Reco", required=False, help="Input histograms")
  parser.add_argument("--particle", dest='particle', default="electron", required=False, help="Particle (electron or pion)")
  parser.add_argument("--plots",dest='plots', nargs='+', type=int, required=True, help='Options to draw different plots')
  args = parser.parse_args()
  plots = args.plots
  plotToString = [
   'Draw all plots',                                                               # 0
   'Draw HCal response(resolution) vs energy plots',                               # 1
   'Draw HCal response(resolution) vs theta plots',                                # 2
   'Draw ECal+HCal TopoClusters response(resolution) vs energy plots',             # 3
   'Draw ECal+HCal TopoClusters response(resolution) vs theta plots',              # 4
  ]
  for i in range(len(plotToString)):
    print ('%s:  ' % i, plotToString[i] )

  print('Now we will:')
  for i in range(len(plots)):
    print ( plotToString[plots[i]] )

  particle = {'electron':'e^{-}','pion':'#pi^{-}','K0L':'K^{0}_{L}'}
  ########################################

  #######################################################################################
  ##  HCal-standalone response linearity and resolution plots @ theta = 68 (38) degree ##
  #######################################################################################
  if 1 in plots or 0 in plots:
    for theta in [68,38]:
      graphs = []
      graphs_sigma = []
      for p in ['electron','pion','K0L']:
        x  = np.array([2, 3, 4, 5, 7, 10, 13, 20, 30, 50, 80, 100, 150, 180],dtype='float64')
        ex = np.zeros(len(x),dtype='float64')  # x-errors
        mu = np.zeros(len(x),dtype='float64')
        emu = np.zeros(len(x),dtype='float64')  # y-errors
        sigma  = np.zeros(len(x),dtype='float64')
        esigma = np.zeros(len(x),dtype='float64')
        for energy,bin in zip(x,range(len(x))):
          # NOTE: electrons are not simulated for Endcap
          #if theta==38 and p=='electron': continue
          inputFile = root.TFile.Open('root://128.141.173.81:1094//home/data/FCC/Reco/HighStat/ALLEGRO_HCalOnly_k4run_%s_E%d_theta%s.root'%(p,energy,theta),"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw('Sum$(CaloTopoClusters.energy)','','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])

          hist = root.TH1F(f'hist_{energy}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw('Sum$(CaloTopoClusters.energy)>>'+hist.GetName())
          mu[bin] = hist.GetMean()
          emu[bin] = hist.GetMeanError()
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigma[bin] = r.n
          esigma[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['HCal-standalone, B = 0T',f'{energy} GeV {particle[p]} @ #theta={theta}#circ','Topo clusters @ EM-scale','','Mean = %.2f GeV' % hist.GetMean()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = '%s_HCalOnly_theta%d_Ebins.pdf('%(p,theta)
          if bin == len(x)-1: opt.Name = '%s_HCalOnly_theta%d_Ebins.pdf)'%(p,theta)
          Helpers.plot(hist,opt)
          mu[bin] = mu[bin]/energy
          emu[bin] = emu[bin]/energy
          # end of energy bin loop
        gr = root.TGraphErrors(len(x), x, mu, ex, emu)
        gr.SetTitle(particle[p])
        graphs.append(gr)
        gr = root.TGraphErrors(len(x), x, sigma, ex, esigma)
        gr.SetTitle(particle[p])
        graphs_sigma.append(gr)
        # end of particles loop

      opt = option()
      opt.xTitle = 'E^{Gen} [GeV]'
      opt.yTitle = '#LTE^{Reco}#GT/E^{Gen}'
      opt.Logx = True
      opt.xRange = [1,200]
      opt.yRange = [0.749,1.2]
      opt.Labels =  ['HCal-standalone, B = 0T',f'#theta={theta}#circ','Topo clusters @ EM-scale']
      opt.Name = f'Linearity_HCalOnly_theta{theta}_Ebins.pdf'
      #if theta==38: Helpers.plotGraph(graphs[1:],opt)
      #else: Helpers.plotGraph(graphs,opt)
      Helpers.plotGraph(graphs,opt)

      func = root.TF1('func_pi','sqrt( ([0]/sqrt(x))*([0]/sqrt(x)) + [1]*[1] )',0.5, 200)
      func.SetLineColor(Helpers.Color(0))
      graphs_sigma[1].Fit(func,'','',0.5, 200)
      opt = option()
      opt.xTitle = 'E^{Gen} [GeV]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.Logx = True
      opt.xRange = [1,200]
      opt.yRange = [0.,0.6]
      opt.Labels =  ['HCal-standalone, B = 0T','#pi^{-} @ #theta=%s#circ'%theta,'Topo clusters @ EM-scale']
      opt.Function = [func,'#frac{%.1f%% }{#sqrt{E}} #oplus %.1f%% ' % (func.GetParameter(0)*100,func.GetParameter(1)*100)]
      opt.Name = f'Resolution_HCalOnly_theta{theta}_pion.pdf'
      Helpers.plotGraph(graphs_sigma[1:2],opt)

      func = root.TF1('func_k','sqrt( ([0]/sqrt(x))*([0]/sqrt(x)) + [1]*[1] )',0.5, 200)
      func.SetLineColor(Helpers.Color(0))
      graphs_sigma[2].Fit(func,'','',0.5, 200)
      opt = option()
      opt.xTitle = 'E^{Gen} [GeV]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.Logx = True
      opt.xRange = [1,200]
      opt.yRange = [0.,0.6]
      opt.Labels =  ['HCal-standalone, B = 0T','K^{0}_{L} @ #theta=%s#circ'%theta,'Topo clusters @ EM-scale']
      opt.Function = [func,'#frac{%.1f%% }{#sqrt{E}} #oplus %.1f%% ' % (func.GetParameter(0)*100,func.GetParameter(1)*100)]
      opt.Name = f'Resolution_HCalOnly_theta{theta}_kaon.pdf'
      Helpers.plotGraph(graphs_sigma[2:3],opt)


  ########################################
  ##  HCal-standalone response vs theta ##
  ########################################
  if 2 in plots or 0 in plots:
    graphs = []
    graphs_sigma = []
    for p in ['electron','pion','K0L']:
      x  = np.array([88, 83, 78, 73, 68, 63, 58, 53, 48, 43, 38, 33, 28],dtype='float64')
      ex = np.zeros(len(x),dtype='float64')  # x-errors
      mu = np.zeros(len(x),dtype='float64')
      emu = np.zeros(len(x),dtype='float64')  # y-errors
      sigma  = np.zeros(len(x),dtype='float64')
      esigma = np.zeros(len(x),dtype='float64')
      for theta,bin in zip(x,range(len(x))):
        inputFile = root.TFile.Open('root://128.141.173.81:1094//home/data/FCC/Reco/HighStat/ALLEGRO_HCalOnly_k4run_%s_E100_theta%d.root'%(p,theta),"read")
        tree = inputFile.Get('events')
        tree.SetEstimate(-1)
        N = tree.Draw('Sum$(CaloTopoClusters.energy)','','goff')
#        N = tree.Draw('(Sum$(HCalEndcapReadout.energy)+Sum$(HCalBarrelReadout.energy))*31.19','','goff')

        values = []
        for i in range(N):
          values.append(tree.GetV1()[i])

        hist = root.TH1F(f'hist_{theta}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
        tree.Draw('Sum$(CaloTopoClusters.energy)>>'+hist.GetName())
#        tree.Draw('(Sum$(HCalEndcapReadout.energy)+Sum$(HCalBarrelReadout.energy))*31.19>>'+hist.GetName())
        mu[bin] = hist.GetMean()
        emu[bin] = hist.GetMeanError()
        hist.Scale(1./hist.Integral(0,-1))
        m = [0.,0.]
        s = [0.,0.]
        Helpers.FitGaus(hist, [1.5,2.], m, s)
        #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
        func = hist.GetFunction('gaus')
        r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
        sigma[bin] = r.n
        esigma[bin] = r.s

        opt = option()
        opt.xTitle = 'E^{reco} [GeV]'
        opt.yTitle = 'Normalized Entries'
        opt.Labels =  ['HCal-standalone',f'100 GeV {particle[p]} @ #theta={theta}#circ','Topo clusters @ EM-scale','','Mean = %.2f GeV' % hist.GetMean()]
        opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
        opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
        opt.Function = hist.GetFunction('gaus')
        opt.Name = 'HCalOnly_%s_E100.pdf('%(p)
        if bin == len(x)-1: opt.Name = 'HCalOnly_%s_E100.pdf)'%(p)
        Helpers.plot(hist,opt)
        mu[bin] = mu[bin]/100.
        emu[bin] = emu[bin]/100.
        # end of energy bin loop
      gr = root.TGraphErrors(len(x), x, mu, ex, emu)
      gr.SetTitle(particle[p])
      graphs.append(gr)
      gr = root.TGraphErrors(len(x), x, sigma, ex, esigma)
      gr.SetTitle(particle[p])
      graphs_sigma.append(gr)
      # end of particles loop

    opt = option()
    opt.xTitle = '#theta [#circ]'
    opt.yTitle = '#LTE^{Reco}#GT/E^{Gen}'
    opt.xRange = [25,90]
    opt.yRange = [0.7,1.3]
    opt.Labels =  ['HCal-standalone','100 GeV','Topo clusters @ EM-scale']
    opt.Name = 'Linearity_vs_theta_HCalOnly.pdf'
    Helpers.plotGraph(graphs,opt)
    opt = option()
    opt.xTitle = '#theta [#circ]'
    opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
    opt.xRange = [25,90]
    opt.yRange = [0.,0.15]
    opt.Labels =  ['HCal-standalone','#pi^{-} @ 100 GeV','Topo clusters @ EM-scale']
    opt.Name = 'Resolution_vs_theta_HCalOnly_pion.pdf'
    Helpers.plotGraph(graphs_sigma[1:2],opt)
    opt = option()
    opt.xTitle = '#theta [#circ]'
    opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
    opt.xRange = [25,90]
    opt.yRange = [0.,0.15]
    opt.Labels =  ['HCal-standalone','K^{0}_{L} @ 100 GeV','Topo clusters @ EM-scale']
    opt.Name = 'Resolution_vs_theta_HCalOnly_kaon.pdf'
    Helpers.plotGraph(graphs_sigma[2:3],opt)


  ##################################################
  ##  ECal+HCal topo-calusters response vs energy ##
  ##################################################
  if 3 in plots or 0 in plots:
    for theta in [68,33]:
      graphs = []
      graphs_sigma = []
      for p in ['pion','K0L']:
        x  = np.array([2, 3, 4, 5, 7, 10, 13, 20, 30, 50, 80, 100, 150, 180],dtype='float64')
        ex = np.zeros(len(x),dtype='float64')  # x-errors
        mu = np.zeros(len(x),dtype='float64')
        emu = np.zeros(len(x),dtype='float64')  # y-errors
        sigma  = np.zeros(len(x),dtype='float64')
        esigma = np.zeros(len(x),dtype='float64')
        for energy,bin in zip(x,range(len(x))):
          inputFile = root.TFile.Open('root://128.141.173.81:1094//home/data/FCC/Reco/HighStat/ALLEGRO_k4run_%s_E%d_theta%s.root'%(p,energy,theta),"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw('Sum$(CaloTopoClusters.energy)','','goff')
#          N = tree.Draw('Sum$(ECalEndcapTurbinePositioned.energy)+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)','','goff')
#          N = tree.Draw('Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy)','','goff')
#          N = tree.Draw('Sum$(ECalEndcapTurbine.energy)+Sum$(ECalBarrelModuleThetaMerged.energy)','','goff')
#          N = tree.Draw('Sum$(HCalEndcapReadout.energy)+Sum$(HCalBarrelReadout.energy)','','goff')


          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])

          hist = root.TH1F(f'hist_{energy}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw('Sum$(CaloTopoClusters.energy)>>'+hist.GetName())
#          tree.Draw('Sum$(ECalEndcapTurbinePositioned.energy)+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)>>'+hist.GetName())
#          tree.Draw('Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy)>>'+hist.GetName())
#          tree.Draw('Sum$(ECalEndcapTurbine.energy)+Sum$(ECalBarrelModuleThetaMerged.energy)>>'+hist.GetName())
#          tree.Draw('Sum$(HCalEndcapReadout.energy)+Sum$(HCalBarrelReadout.energy)>>'+hist.GetName())

          mu[bin] = hist.GetMean()
          emu[bin] = hist.GetMeanError()
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,2.], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigma[bin] = r.n
          esigma[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal','Topo clusters @ EM-scale', f'{energy} GeV {particle[p]} @ #theta=68#circ','','Mean = %.2f GeV' % hist.GetMean()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = '%s_theta%s.pdf('%(p,theta)
          if bin == len(x)-1: opt.Name = '%s_theta%s.pdf)'%(p,theta)
          Helpers.plot(hist,opt)
          mu[bin] = mu[bin]/energy
          emu[bin] = emu[bin]/energy
          # end of energy bin loop
        gr = root.TGraphErrors(len(x), x, mu, ex, emu)
        gr.SetTitle(particle[p])
        graphs.append(gr)
        gr = root.TGraphErrors(len(x), x, sigma, ex, esigma)
        gr.SetTitle(particle[p])
        graphs_sigma.append(gr)
        # end of particles loop

      opt = option()
      opt.xTitle = 'E^{Gen} [GeV]'
      opt.yTitle = '#LTE^{Reco}#GT/E^{Gen}'
      opt.Logx = True
      opt.xRange = [1,300]
      opt.yRange = [0.68,1.25]
#      opt.yRange = [0.,1.2]
      opt.Labels =  ['ECal+HCal',f'#theta={theta}#circ','Topo clusters @ EM-scale']
#      opt.Labels =  ['ECal+HCal',f'#theta={theta}#circ','HCal hits']
      opt.Name = f'Linearity_theta{theta}.pdf'
      Helpers.plotGraph(graphs,opt)

      func = root.TF1('func_pi','sqrt( ([0]/sqrt(x))*([0]/sqrt(x)) + [1]*[1] )',0.5, 200)
      func.SetLineColor(Helpers.Color(0))
      graphs_sigma[0].Fit(func,'','',0.5, 200)
      opt = option()
      opt.xTitle = 'E^{Gen} [GeV]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.Logx = True
      opt.xRange = [1,200]
      opt.yRange = [0.,0.6]
      opt.Labels =  ['ECal+HCal','#pi^{-} @ #theta=%s#circ'%theta,'Topo clusters @ EM-scale']
      opt.Function = [func,'#frac{%.1f%% }{#sqrt{E}} #oplus %.1f%% ' % (func.GetParameter(0)*100,func.GetParameter(1)*100)]
      opt.Name = f'Resolution_pion_theta{theta}.pdf'
      Helpers.plotGraph(graphs_sigma[0:1],opt)

      func = root.TF1('func_k','sqrt( ([0]/sqrt(x))*([0]/sqrt(x)) + [1]*[1] )',0.5, 200)
      func.SetLineColor(Helpers.Color(0))
      graphs_sigma[1].Fit(func,'','',0.5, 200)
      opt = option()
      opt.xTitle = 'E^{Gen} [GeV]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.Logx = True
      opt.xRange = [1,200]
      opt.yRange = [0.,0.6]
      opt.Labels =  ['ECal+HCal','K^{0}_{L} @ #theta=%s#circ'%theta,'Topo clusters @ EM-scale']
      opt.Function = [func,'#frac{%.1f%% }{#sqrt{E}} #oplus %.1f%% ' % (func.GetParameter(0)*100,func.GetParameter(1)*100)]
      opt.Name = f'Resolution_kaon_theta{theta}.pdf'
      Helpers.plotGraph(graphs_sigma[1:2],opt)


  ##################################
  ##  ECal+HCal response vs theta ##
  ##################################
  if 4 in plots or 0 in plots:
    graphs = []
    graphs_sigma = []
    for p in ['pion','K0L']:
      x  = np.array([88, 83, 78, 73, 68, 63, 58, 53, 48, 43, 38, 33, 28],dtype='float64')
#      x  = np.array([88, 83, 78, 73, 68, 63, 58],dtype='float64')

      ex = np.zeros(len(x),dtype='float64')  # x-errors
      mu = np.zeros(len(x),dtype='float64')
      emu = np.zeros(len(x),dtype='float64')  # y-errors
      sigma  = np.zeros(len(x),dtype='float64')
      esigma = np.zeros(len(x),dtype='float64')
      for theta,bin in zip(x,range(len(x))):
        inputFile = root.TFile.Open('root://128.141.173.81:1094//home/data/FCC/Reco/HighStat/ALLEGRO_k4run_%s_E100_theta%d.root'%(p,theta),"read")
        tree = inputFile.Get('events')
        tree.SetEstimate(-1)
        N = tree.Draw('Sum$(CaloTopoClusters.energy)','','goff')
#        N = tree.Draw('Sum$(ECalEndcapTurbinePositioned.energy)+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)','','goff')
#        N = tree.Draw('Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy)','','goff')
#        N = tree.Draw('Sum$(ECalEndcapTurbine.energy)+Sum$(ECalBarrelModuleThetaMerged.energy)','','goff')
#        N = tree.Draw('Sum$(HCalEndcapReadout.energy)+Sum$(HCalBarrelReadout.energy)','','goff')

        values = []
        for i in range(N):
          values.append(tree.GetV1()[i])

        hist = root.TH1F(f'hist_{theta}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
        tree.Draw('Sum$(CaloTopoClusters.energy)>>'+hist.GetName())
#        tree.Draw('Sum$(ECalEndcapTurbinePositioned.energy)+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)>>'+hist.GetName())
#        tree.Draw('Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy)>>'+hist.GetName())
#        tree.Draw('Sum$(ECalEndcapTurbine.energy)+Sum$(ECalBarrelModuleThetaMerged.energy)>>'+hist.GetName())
#        tree.Draw('Sum$(HCalEndcapReadout.energy)+Sum$(HCalBarrelReadout.energy)>>'+hist.GetName())


        mu[bin] = hist.GetMean()
        emu[bin] = hist.GetMeanError()
        hist.Scale(1./hist.Integral(0,-1))
        m = [0.,0.]
        s = [0.,0.]
        Helpers.FitGaus(hist, [1.5,2.], m, s)
        #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
        func = hist.GetFunction('gaus')
        r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
        sigma[bin] = r.n
        esigma[bin] = r.s

        opt = option()
        opt.xTitle = 'E^{reco} [GeV]'
        opt.yTitle = 'Normalized Entries'
        opt.Labels =  ['ECal+HCal',f'100 GeV {particle[p]} @ #theta={theta}#circ','Topo clusters @ EM-scale','','Mean = %.2f GeV' % hist.GetMean()]
        opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
        opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
#        opt.Function = hist.GetFunction('gaus')
        opt.Name = '%s_E100.pdf('%(p)
        if bin == len(x)-1: opt.Name = '%s_E100.pdf)'%(p)
        Helpers.plot(hist,opt)
        mu[bin] = mu[bin]/100.
        emu[bin] = emu[bin]/100.
        # end of energy bin loop
      gr = root.TGraphErrors(len(x), x, mu, ex, emu)
      gr.SetTitle(particle[p])
      graphs.append(gr)
      gr = root.TGraphErrors(len(x), x, sigma, ex, esigma)
      gr.SetTitle(particle[p])
      graphs_sigma.append(gr)
      # end of particles loop

    opt = option()
    opt.xTitle = '#theta [#circ]'
    opt.yTitle = '#LTE^{Reco}#GT/E^{Gen}'
    opt.xRange = [25,90]
#    opt.yRange = [0.7,1.3]
    opt.yRange = [0.0,5.3]
    opt.Labels =  ['ECal+HCal','100 GeV','Topo clusters @ EM-scale']
#    opt.Labels =  ['ECal+HCal','100 GeV','ECal cells @ EM-scale']
#    opt.Labels =  ['ECal+HCal','100 GeV','HCal cells @ EM-scale']
#    opt.Labels =  ['ECal+HCal','100 GeV','HCal hits']
    opt.Name = 'Linearity_vs_theta.pdf'
    Helpers.plotGraph(graphs,opt)

    opt = option()
    opt.xTitle = '#theta [#circ]'
    opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
    opt.xRange = [25,90]
    opt.yRange = [0.,0.15]
    opt.Labels =  ['ECal+HCal','#pi^{-} @ 100 GeV','Topo clusters @ EM-scale']
    opt.Name = 'Resolution_vs_theta_pion.pdf'
    Helpers.plotGraph(graphs_sigma[0:1],opt)
    opt = option()
    opt.xTitle = '#theta [#circ]'
    opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
    opt.xRange = [25,90]
    opt.yRange = [0.,0.15]
    opt.Labels =  ['HCal-standalone','K^{0}_{L} @ 100 GeV','Topo clusters @ EM-scale']
    opt.Name = 'Resolution_vs_theta_kaon.pdf'
    Helpers.plotGraph(graphs_sigma[1:2],opt)
