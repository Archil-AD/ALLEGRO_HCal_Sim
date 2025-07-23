# ALLEGRO_HCal_Sim
## checkout
```
git clone https://github.com/Archil-AD/ALLEGRO_HCal_Sim.git
cd ALLEGRO_HCal_Sim/
```
## FCC software setup 
```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```
## Run HCal-standalone simulation and reconstruction
### setup
Only for the first time do the following:
```
git clone https://github.com/key4hep/k4geo.git
cd k4geo/
export K4GEO=$PWD
cp ../run/ALLEGRO_o1_v03_HCal_standalone.xml $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/
cd ../run
```
Next time you need to do only the following:
```
cd k4geo/
export K4GEO=$PWD
cd ../run
```
### simulation
```
ddsim --enableGun --gun.distribution uniform --gun.energy "100*GeV" --gun.thetaMin "68*deg" --gun.thetaMax "68*deg" --gun.particle e- --numberOfEvents 500 --outputFile ALLEGRO_Sim_HCalOnly_E100_theta68.root --random.enableEventSeed --random.seed 4352 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03_HCal_standalone.xml
```
### reconstruction
```
k4run run_reco_HCal.py --inputFiles ALLEGRO_Sim_HCalOnly_E100_theta68.root --outputFile ALLEGRO_Rec_HCalOnly_E100_theta68.root
```
