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
  parser.add_argument("--inputDir", dest='inputDir', nargs='+', default=[""], required=True, help="Input histograms")
  parser.add_argument("--outputDir", dest='outputDir', default="", required=False, help="Output plots directory")
  parser.add_argument("--energy",dest='energy', nargs='+', default=[30, 10], type=int, required=False, help='Energy')
  parser.add_argument("--theta",dest='theta', nargs='+', default=[88, 78, 68, 58], type=int, required=False, help='Theta')
  parser.add_argument("--alpha",dest='alpha', nargs='+', default=[5, 10, 15, 25], type=int, required=False, help='Opening angle alpha')
  parser.add_argument("--plot",dest='plot', nargs='+', type=int, default=[0], required=False, help='Options to draw different plots')
  args = parser.parse_args()
  plots = args.plot
  plotToString = [
   'Draw all plots',                                                               # 0
   'Draw calo energy res vs theta',                                                # 1
   'Draw PFO energy res vs theta',                                                 # 2
   'Draw Perfect PFO energy res vs theta',                                         # 3
   'Draw energy res vs alpha @ theta=68',                                          # 4
   'Draw energy res vs alpha @ theta=48',                                          # 5
   'Draw PFO energy comparison',                                                   # 6
  ]
  for i in range(len(plotToString)):
    print ('%s:  ' % i, plotToString[i] )

  print('Now we will:')
  for i in range(len(plots)):
    print ( plotToString[plots[i]] )

  # create output directory inside 'plots'
  if not os.path.exists(f'plots/{args.outputDir}'):
    os.makedirs(f'plots/{args.outputDir}')


  #######################################################
  ##  Calorimeter resolution for two particles events  ##
  #######################################################
  if 1 in plots or 0 in plots:
      graphs_sigma = []
      graphs_gauss_sigma = []
      caloEnergy = '(Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy))*1.15+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)'
      alpha = np.array(args.alpha,dtype='float64')
      for a in alpha:
        theta = np.array(args.theta,dtype='float64')
        sigma  = np.zeros(len(theta),dtype='float64')
        esigma = np.zeros(len(theta),dtype='float64')
        sigmaFit  = np.zeros(len(theta),dtype='float64')
        esigmaFit = np.zeros(len(theta),dtype='float64')
        for th,bin in zip(args.theta,range(len(theta))):
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(caloEnergy,'','goff')

          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(caloEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{calo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]

          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/CaloEnergy_alpha{a:.1f}.pdf('
          if bin == len(theta)-1: opt.Name = f'{args.outputDir}/CaloEnergy_alpha{a:.1f}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(theta), theta, sigma, np.zeros(len(theta),dtype='float64'), esigma)
        gr.SetTitle(f'#alpha={a:.0f}#circ')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(theta), theta, sigmaFit, np.zeros(len(theta),dtype='float64'), esigmaFit)
        gr.SetTitle(f'#alpha={a:.0f}#circ')
        graphs_gauss_sigma.append(gr)

      opt = option()
      opt.xTitle = '#theta^{#pi^{-}} [#circ]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.xRange = [55,90]
      opt.yRange = [0.01,0.16]
      opt.Labels = ['ECal+HCal','30 GeV #pi^{-} + 10 GeV K^{0}_{L}']
      opt.Labels += ['Calorimeter energy']
      opt.Name = f'{args.outputDir}/Calo_resolution_vs_theta.pdf'
      Helpers.plotGraph(graphs_gauss_sigma,opt)

  #######################################################
  ##      PFO resolution for two particles events      ##
  #######################################################
  if 2 in plots or 0 in plots:
      graphs_sigma = []
      graphs_gauss_sigma = []
      pfoEnergy = 'Sum$(PandoraPFANewPFOs.energy)'
      alpha = np.array(args.alpha,dtype='float64')
      for a in alpha:
        theta = np.array(args.theta,dtype='float64')
        sigma  = np.zeros(len(theta),dtype='float64')
        esigma = np.zeros(len(theta),dtype='float64')
        sigmaFit  = np.zeros(len(theta),dtype='float64')
        esigmaFit = np.zeros(len(theta),dtype='float64')
        for th,bin in zip(args.theta,range(len(theta))):
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy,'','goff')

          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]

          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/PFOEnergy_alpha{a:.1f}.pdf('
          if bin == len(theta)-1: opt.Name = f'{args.outputDir}/PFOEnergy_alpha{a:.1f}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(theta), theta, sigma, np.zeros(len(theta),dtype='float64'), esigma)
        gr.SetTitle(f'#alpha={a:.1f}#circ')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(theta), theta, sigmaFit, np.zeros(len(theta),dtype='float64'), esigmaFit)
        gr.SetTitle(f'#alpha={a:.1f}#circ')
        graphs_gauss_sigma.append(gr)

      opt = option()
      opt.xTitle = '#theta^{#pi^{-}} [#circ]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.xRange = [55,90]
      opt.yRange = [0.01,0.16]
      opt.Labels = ['ECal+HCal','30 GeV #pi^{-} + 10 GeV K^{0}_{L}']
      opt.Labels += ['PFO energy']
      opt.Name = f'{args.outputDir}/PFO_resolution_vs_theta.pdf'
      Helpers.plotGraph(graphs_gauss_sigma,opt)

  #######################################################
  ##  Pefect PFO resolution for two particles events   ##
  #######################################################
  if 3 in plots or 0 in plots:
      graphs_sigma = []
      graphs_gauss_sigma = []
      pfoEnergy = 'Sum$(PerfectPFOs.energy)'
      alpha = np.array(args.alpha,dtype='float64')
      for a in [5.0]:
        theta = np.array(args.theta,dtype='float64')
        sigma  = np.zeros(len(theta),dtype='float64')
        esigma = np.zeros(len(theta),dtype='float64')
        sigmaFit  = np.zeros(len(theta),dtype='float64')
        esigmaFit = np.zeros(len(theta),dtype='float64')
        for th,bin in zip(args.theta,range(len(theta))):
          if not os.path.exists(args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy,'','goff')

          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]

          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/PerfectPFOEnergy_alpha{a:.1f}.pdf('
          if bin == len(theta)-1: opt.Name = f'{args.outputDir}/PerfectPFOEnergy_alpha{a:.1f}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(theta), theta, sigma, np.zeros(len(theta),dtype='float64'), esigma)
        gr.SetTitle(f'#alpha={a:.1f}#circ')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(theta), theta, sigmaFit, np.zeros(len(theta),dtype='float64'), esigmaFit)
        gr.SetTitle(f'#alpha={a:.1f}#circ')
        graphs_gauss_sigma.append(gr)

      opt = option()
      opt.xTitle = '#theta^{#pi^{-}} [#circ]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.xRange = [55,90]
      opt.yRange = [0.01,0.16]
      opt.Labels = ['ECal+HCal','30 GeV #pi^{-} + 10 GeV K^{0}_{L}','#alpha=5#circ']
      opt.Labels += ['Perfect PFO energy']
      opt.Name = f'{args.outputDir}/PerfectPFO_resolution_vs_theta.pdf'
      Helpers.plotGraph(graphs_gauss_sigma,opt)


  #######################################################
  ##     Energy resolution vs alpha for theta = 68     ##
  #######################################################
  if 4 in plots or 0 in plots:
      graphs_sigma = []
      graphs_gauss_sigma = []
      alpha = np.array(args.alpha,dtype='float64')
      for th in [68]:
        ######## Calo
        caloHists = []
        caloEnergy = '(Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy))*1.15+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)'
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'calo_hist_{th}_alpha{a}','Calo',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(caloEnergy+'>>'+hist.GetName(),'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          caloHists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(caloEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/CaloEnergy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/CaloEnergy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'Calo')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'Calo')
        graphs_gauss_sigma.append(gr)

        ######## PFO
        pfoHists = []
        pfoEnergy = 'Sum$(PandoraPFANewPFOs.energy)'
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'pfo_hist_{th}_alpha{a}','PFO',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy+'>>'+hist.GetName(),'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          pfoHists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/PFOEnergy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/PFOEnergy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'PFO')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'PFO')
        graphs_gauss_sigma.append(gr)

        ######## perfect PFO
        perfectPfoHists = []
        pfoEnergy = 'Sum$(PerfectPFOs.energy)'
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'perfectPfo_hist_{th}_alpha{a}','Perfect PFO',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'_perfectPFA/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy,'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          perfectPfoHists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/PerfectPFOEnergy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/PerfectPFOEnergy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'Perfect PFO')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'Perfect PFO')
        graphs_gauss_sigma.append(gr)

        # plot superimposed distributions for Calo, PFO and perfect PFO
        for a,i in zip(args.alpha,range(len(alpha))):
          opt = option()
          opt.xTitle = 'Reco energy [GeV]'
          opt.yTitle = 'Events'
          #opt.xRange = [20,70]
          #opt.yRange = [0.01,0.26]
          opt.Labels = ['ECal+HCal','30 GeV #pi^{-}@%d#circ + 10 GeV K^{0}_{L}'%th,f'#alpha={a}#circ']
          opt.Name = f'{args.outputDir}/Hist_theta{th}_alpha{a}.pdf'
          opt.Hists = [caloHists[i],pfoHists[i],perfectPfoHists[i]]
          dummyHist = root.TH1F(f'dummy_hist_{th}_{a}','',50,20,70)
          dummyHist.SetMaximum(perfectPfoHists[i].GetMaximum())
          Helpers.plot(dummyHist,opt)

      opt = option()
      opt.xTitle = '#alpha [#circ]'
      opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
      opt.xRange = [0,30]
      opt.yRange = [0.01,0.26]
      opt.Labels = ['ECal+HCal','30 GeV #pi^{-} + 10 GeV K^{0}_{L}',f'#theta={th}#circ']
      opt.Name = f'{args.outputDir}/Resolution_vs_alpha.pdf'
      Helpers.plotGraph(graphs_gauss_sigma,opt)


  #######################################################
  ##     Energy resolution vs alpha for theta = 48     ##
  #######################################################
  if 5 in plots or 0 in plots:
      alpha = np.array(args.alpha,dtype='float64')
      for th in [48]:
        graphs_sigma = []
        graphs_gauss_sigma = []
        ######## Calo
        caloHists = []
        caloEnergy = '(Sum$(HCalEndcapReadoutPositioned.energy)+Sum$(HCalBarrelReadoutPositioned.energy))*1.15+Sum$(ECalBarrelModuleThetaMergedPositioned.energy)'
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'calo_hist_{th}_alpha{a}','Calo',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(caloEnergy+'>>'+hist.GetName(),'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          caloHists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(caloEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/CaloEnergy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/CaloEnergy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'Calo')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'Calo')
        graphs_gauss_sigma.append(gr)

        ######## PFO
        pfoHists = []
        pfoEnergy = 'Sum$(PandoraPFANewPFOs.energy)'
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'pfo_hist_{th}_alpha{a}','PFO',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy+'>>'+hist.GetName(),'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          pfoHists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/PFOEnergy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/PFOEnergy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'PFO')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'PFO')
        graphs_gauss_sigma.append(gr)

        # plot superimposed distributions for Calo, PFO and perfect PFO
        for a,i in zip(args.alpha,range(len(alpha))):
          opt = option()
          opt.xTitle = 'Reco energy [GeV]'
          opt.yTitle = 'Events'
          #opt.xRange = [20,70]
          #opt.yRange = [0.01,0.26]
          opt.Labels = ['ECal+HCal','30 GeV #pi^{-}@%d#circ + 10 GeV K^{0}_{L}'%th,f'#alpha={a}#circ']
          opt.Name = f'{args.outputDir}/Hist_theta{th}_alpha{a}.pdf'
          opt.Hists = [caloHists[i],pfoHists[i]]
          dummyHist = root.TH1F(f'dummy_hist_{th}_{a}','',50,20,70)
          dummyHist.SetMaximum(max([pfoHists[i].GetMaximum(),caloHists[i].GetMaximum()]))
          Helpers.plot(dummyHist,opt)

        opt = option()
        opt.xTitle = '#alpha [#circ]'
        opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
        opt.xRange = [0,30]
        opt.yRange = [0.01,0.26]
        opt.Labels = ['ECal+HCal','30 GeV #pi^{-} + 10 GeV K^{0}_{L}',f'#theta={th}#circ']
        opt.Name = f'{args.outputDir}/Resolution_vs_alpha_theta{th}.pdf'
        Helpers.plotGraph(graphs_gauss_sigma,opt)

  #######################################################
  ##     Energy resolution vs alpha for theta = 48     ##
  #######################################################
  if 6 in plots or 0 in plots:
      if len(args.inputDir)!=2:
        print('Please provide two input directories!')
        sys.exit()

      alpha = np.array(args.alpha,dtype='float64')
      for th in args.theta:
        graphs_sigma = []
        graphs_gauss_sigma = []
        pfoEnergy = 'Sum$(PandoraPFANewPFOs.energy)'
        #pfoEnergy = 'Sum$(PerfectPFOs.energy)'

        ######## pfo0
        pfo0Hists = []
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'pfo0_hist_{th}_alpha{a}','highest granularity',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[0]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy+'>>'+hist.GetName(),'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          pfo0Hists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/pfo0Energy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/pfo0Energy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'pfo0')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'pfo0')
        graphs_gauss_sigma.append(gr)

        ######## pfo1
        pfo1Hists = []
        alpha = np.array(args.alpha,dtype='float64')
        sigma  = np.zeros(len(alpha),dtype='float64')
        esigma = np.zeros(len(alpha),dtype='float64')
        sigmaFit  = np.zeros(len(alpha),dtype='float64')
        esigmaFit = np.zeros(len(alpha),dtype='float64')
        for a,bin in zip(args.alpha,range(len(alpha))):
          # histogram with fixed range and bins for comparison between Calo, PFO and perfect PFO
          hist = root.TH1F(f'pfo0_hist_{th}_alpha{a}','Coarser granularity',50,20,70)
          hist.SetDirectory(0)
          if not os.path.exists(args.inputDir[1]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root'):
            print('File does not exist',args.inputDir[1]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
            continue
          print('Opening:',args.inputDir[1]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root')
          inputFile = root.TFile.Open(args.inputDir[1]+f'/ALLEGRO_k4run_alpha{a:.1f}_K0L_E30_pion_E10_theta{th}.root',"read")
          tree = inputFile.Get('events')
          tree.SetEstimate(-1)
          N = tree.Draw(pfoEnergy+'>>'+hist.GetName(),'','goff')
          values = []
          for i in range(N):
            values.append(tree.GetV1()[i])
            hist.Fill(tree.GetV1()[i])
          pfo1Hists.append(hist)

          hist = root.TH1F(f'hist_{th}','',80,np.mean(values) - 4*np.std(values), np.mean(values) + 6*np.std(values))
          tree.Draw(pfoEnergy+'>>'+hist.GetName())
          r = unc.ufloat(hist.GetRMS(),hist.GetRMSError())/unc.ufloat(hist.GetMean(),hist.GetMeanError())
          sigma[bin] = r.n
          esigma[bin] = r.s
          hist.Scale(1./hist.Integral(0,-1))
          m = [0.,0.]
          s = [0.,0.]
          Helpers.FitGaus(hist, [1.5,1.5], m, s)
          #hist.Fit('gaus','QIM','',np.mean(values) - 1.5*np.std(values), np.mean(values) + 2*np.std(values))
          func = hist.GetFunction('gaus')
          r = unc.ufloat(func.GetParameter(2),func.GetParError(2))/unc.ufloat(func.GetParameter(1),func.GetParError(1))
          sigmaFit[bin] = r.n
          esigmaFit[bin] = r.s

          opt = option()
          opt.xTitle = 'E^{reco}_{pfo} [GeV]'
          opt.yTitle = 'Normalized Entries'
          opt.Labels =  ['ECal+HCal',f'30 GeV #pi^{{-}}@#theta={th}#circ + 10 GeV K^{{0}}_{{L}}',f'#alpha={a:.0f}#circ']
          opt.Labels += ['Mean = %.2f GeV' % hist.GetMean(),'Std.Dev. = %.2f GeV' % hist.GetRMS()]
          opt.Labels += ['#mu^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(1),func.GetParError(1))]
          opt.Labels += ['#sigma^{Gauss} = %.2f#pm %.2f GeV' % (func.GetParameter(2),func.GetParError(2))]
          opt.Labels += ['#sigma^{Gauss}/#mu^{Gauss} = %s' % r]
          opt.Function = hist.GetFunction('gaus')
          opt.Name = f'{args.outputDir}/pfo1Energy_theta{th}.pdf('
          if bin == len(alpha)-1: opt.Name = f'{args.outputDir}/pfo1Energy_theta{th}.pdf)'
          Helpers.plot(hist,opt)
        gr = root.TGraphErrors(len(alpha), alpha, sigma, np.zeros(len(alpha),dtype='float64'), esigma)
        gr.SetTitle(f'pfo1')
        graphs_sigma.append(gr)
        gr = root.TGraphErrors(len(alpha), alpha, sigmaFit, np.zeros(len(alpha),dtype='float64'), esigmaFit)
        gr.SetTitle(f'pfo1')
        graphs_gauss_sigma.append(gr)

        # plot superimposed distributions
        for a,i in zip(args.alpha,range(len(alpha))):
          opt = option()
          opt.xTitle = 'Total PFO energy [GeV]'
          opt.yTitle = 'Events'
          #opt.xRange = [20,70]
          #opt.yRange = [0.01,0.26]
          opt.Labels = ['ECal+HCal','30 GeV #pi^{-}@%d#circ + 10 GeV K^{0}_{L}'%th,f'#alpha={a}#circ']
          opt.Name = f'{args.outputDir}/Hist_theta{th}_alpha{a}.pdf'
          opt.Hists = [pfo0Hists[i],pfo1Hists[i]]
          dummyHist = root.TH1F(f'dummy_hist_{th}_{a}','',50,20,70)
          dummyHist.SetMaximum(max([pfo0Hists[i].GetMaximum(),pfo1Hists[i].GetMaximum()]))
          Helpers.plot(dummyHist,opt)

        opt = option()
        opt.xTitle = '#alpha [#circ]'
        opt.yTitle = '#sigma(E^{Reco})/#mu(E^{Reco})'
        opt.xRange = [0,30]
        opt.yRange = [0.01,0.26]
        opt.Labels = ['ECal+HCal','30 GeV #pi^{-} + 10 GeV K^{0}_{L}',f'#theta={th}#circ']
        opt.Name = f'{args.outputDir}/PFO_comp_res_vs_alpha_theta{th}.pdf'
        Helpers.plotGraph(graphs_gauss_sigma,opt)
