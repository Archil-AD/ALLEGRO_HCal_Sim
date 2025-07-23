#!/usr/bin/python
import argparse
import time
from datetime import date
import os
import sys
import subprocess

parser = argparse.ArgumentParser(description='Process some parameters')
parser.add_argument("--inputDir", dest='inputDir', default='', required=False, help="Input SIM files directory")
parser.add_argument("--outputDir", dest='outputDir', default='', required=True, help="Output directory")
parser.add_argument("--outputDest", dest='outputDest', default='', required=True, help="Output destination")
parser.add_argument("--energy",dest='energy', nargs='+', type=int, required=True, help='Energy values')
parser.add_argument("--theta",dest='theta', nargs='+', type=int, required=True, help='Theta values')
parser.add_argument("--particle", dest='particle', default='', required=True, help="Particle [e-,pi-,kaon0L]")
parser.add_argument("--numberOfEvents", dest='nevts', type=int, default=2500, required=False, help="Number of events")
parser.add_argument("--run", dest='run', default='', required=True, help="ddsim or k4run")
parser.add_argument("--HCalOnly", dest='hcalonly', action='store_true', default=False, required=False, help="Run HCal-standalone")
parser.add_argument("--jobFlavour", dest='jobFlavour', default="workday", required=False, help="HTCondor JobFlavour")
parser.add_argument("--xrd", dest='xrd', action='store_true', default=True, required=False, help="Copy output using XRootD")

args = parser.parse_args()
timestamp = int(time.time())
today = date.today().strftime("%Y-%m-%d")
energy = [str(i) for i in args.energy]
theta = [str(i) for i in args.theta]
particle = args.particle

if args.run == 'k4run' and args.inputDir == '':
  print('Please provide input direcotry containing SIM files')
  sys.exit()

particleDict = {'pi-':'pion','pi+':'pion','e-':'electron','e+':'positron','kaon0L':'K0L'}
particleDictTex = {'pi-':'#pi^{-}','pi+':'#pi^{+}','e-':'e^{-}','e+':'e^{+}','kaon0L':'K^{0}_{L}'}

hcalonly = '_'
if args.hcalonly: hcalonly = '_HCalOnly_'
scriptName = 'empty.sh'

if args.run == 'ddsim':
  outputName = f'ALLEGRO{hcalonly}ddsim_{particleDict[particle]}_E'
  scriptName = f'run_ddsim{hcalonly}{particleDict[particle]}'
  condor_commands = []
  if len(energy) == 1 and len(theta) > 1:
     condor_commands.append('queue %s' % len(theta))
     scriptName += f'_theta_E{energy[0]}.sh'
     bash_commands = [
      '#!/bin/bash',
      'echo "" > empty.file',
      'theta=\"'+' '.join(theta)+'\"',
      'N=$1',
      '((N++))',
      'theta=`echo $theta | cut -d " " -f $N`',
     ]
     if args.hcalonly: bash_commands+=[
       'ddsim --enableGun --gun.distribution uniform --gun.energy "%s*GeV" --gun.thetaMin "${theta}*deg" --gun.thetaMax "${theta}*deg" --gun.particle %s --numberOfEvents %s --outputFile %s_theta${theta}.root --random.enableEventSeed --random.seed 4352${1} --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03_HCal_standalone.xml' % (energy[0],particle,args.nevts,outputName+energy[0]),
     ]
     else: bash_commands+=[
       'ddsim --enableGun --gun.distribution uniform --gun.energy "%s*GeV" --gun.thetaMin "${theta}*deg" --gun.thetaMax "${theta}*deg" --gun.particle %s --numberOfEvents %s --outputFile %s_theta${theta}.root --random.enableEventSeed --random.seed 4352${1} --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml' % (energy[0],particle,args.nevts,outputName+energy[0]),
     ]

  if len(theta) == 1 and len(energy) > 1:
     condor_commands.append('queue %s' % len(energy))
     scriptName += f'_energy_theta{theta[0]}.sh'
     bash_commands = [
      '#!/bin/bash',
      'echo "" > empty.file',
      'energy=\"'+' '.join(energy)+'\"',
      'N=$1',
      '((N++))',
      'energy=`echo $energy | cut -d " " -f $N`',
     ]
     if args.hcalonly: bash_commands+=[
       'ddsim --enableGun --gun.distribution uniform --gun.energy "${energy}*GeV" --gun.thetaMin "%s*deg" --gun.thetaMax "%s*deg" --gun.particle %s --numberOfEvents %s --outputFile %s${energy}_theta%s.root --random.enableEventSeed --random.seed 4352${1} --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03_HCal_standalone.xml' % (theta[0],theta[0],particle,args.nevts,outputName,theta[0]),
     ]
     else: bash_commands+=[
       'ddsim --enableGun --gun.distribution uniform --gun.energy "${energy}*GeV" --gun.thetaMin "%s*deg" --gun.thetaMax "%s*deg" --gun.particle %s --numberOfEvents %s --outputFile %s${energy}_theta%s.root --random.enableEventSeed --random.seed 4352${1} --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml' % (theta[0],theta[0],particle,args.nevts,outputName,theta[0]),
     ]
  bash_commands+= [f'mkdir {args.outputDir} && mv ALLEGRO*ddsim*root {args.outputDir}/']

if args.run == 'k4run':
  inputName = f'{args.inputDir}/ALLEGRO{hcalonly}ddsim_{particleDict[particle]}_E'
  outputName = f'ALLEGRO{hcalonly}k4run_{particleDict[particle]}_E'
  scriptName = f'run_k4run{hcalonly}{particleDict[particle]}'
  condor_commands = []
  if len(energy) == 1 and len(theta) > 1:
     condor_commands.append('queue %s' % len(theta))
     scriptName += f'_theta_E{energy[0]}.sh'
     bash_commands = [
      '#!/bin/bash',
      'echo "" > empty.file',
      'cp -r %s/run ./ && cd run/' % os.getcwd(),
      'theta=\"'+' '.join(theta)+'\"',
      'N=$1',
      '((N++))',
      'theta=`echo $theta | cut -d " " -f $N`',
     ]
     if args.hcalonly: bash_commands+=[
       'k4run run_reco_HCal.py --inputFiles %s_theta${theta}.root --outputFile %s_theta${theta}.root' % (inputName+energy[0],outputName+energy[0]),
     ]
     else: bash_commands+=[
       'k4run run_reco.py --inputFiles %s_theta${theta}.root --outputFile %s_theta${theta}.root' % (inputName+energy[0],outputName+energy[0]),
     ]

  if len(theta) == 1 and len(energy) > 1:
     condor_commands.append('queue %s' % len(energy))
     scriptName += f'_energy_theta{theta[0]}.sh'
     bash_commands = [
      '#!/bin/bash',
      'echo "" > empty.file',
      'cp -r %s/run ./ && cd run/' % os.getcwd(),
      'energy=\"'+' '.join(energy)+'\"',
      'N=$1',
      '((N++))',
      'energy=`echo $energy | cut -d " " -f $N`',
     ]
     if args.hcalonly: bash_commands+=[
       'k4run run_reco_HCal.py --inputFiles %s${energy}_theta%s.root --outputFile %s${energy}_theta%s.root' % (inputName,theta[0],outputName,theta[0]),
     ]
     else: bash_commands+=[
       'k4run run_reco.py --inputFiles %s${energy}_theta%s.root --outputFile %s${energy}_theta%s.root' % (inputName,theta[0],outputName,theta[0]),
     ]
  bash_commands+= [f'mkdir {args.outputDir} && mv ALLEGRO*k4run*root {args.outputDir}/']

if args.xrd: bash_commands+= [f'xrdcp -rf {args.outputDir} {args.outputDest}']
else: bash_commands+= [f'cp -r {args.outputDir} {args.outputDest}']

###################################################
# write to a .sh file
with open(scriptName, "w") as file:
  for cmd in bash_commands:
      file.write(cmd + "\n")
os.chmod(scriptName, 0o755)

condor_commands = [
      f'executable            = {scriptName}',
      'arguments             = $(ProcId)',
      f'output                = output/{scriptName[:-3]}.$(ClusterId).$(ProcId).out',
      f'error                 = error/{scriptName[:-3]}.$(ClusterId).$(ProcId).err',
      f'log                   = log/{scriptName[:-3]}.$(ClusterId).log',
      'transfer_output_files = empty.file',
      'getenv                = True',
      #f'output_destination    = {args.outputDest}',
      f'+JobFlavour           = "{args.jobFlavour}"',
] + condor_commands

# write to a condor submission file
with open(scriptName[:-3]+'.sub', "w") as file:
  for cmd in condor_commands:
      file.write(cmd + "\n")

# submit condor jobs
subprocess.run('condor_submit '+scriptName[:-3]+'.sub', shell=True)
