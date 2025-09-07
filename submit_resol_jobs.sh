## Simulation
if [ "$1" == "ddsim" ]; then
outputDest="root://128.141.173.81:1094//home/data/FCC/Sim/"
outputDir="Resolution"
energy="2 3 4 5 7 10 13 20 30 50 80 100 150 180"
theta="88 83 78 73 63 58 53 48 43 38 33 28 23 18"

# ECal+HCal Barrel
python python/submitJobs.py --energy ${energy} --theta 68 --particle pi- --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow
python python/submitJobs.py --energy ${energy} --theta 68 --particle kaon0L --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow
# HCal Barrel
python python/submitJobs.py --energy ${energy} --theta 68 --particle e- --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --HCalOnly
python python/submitJobs.py --energy ${energy} --theta 68 --particle pi- --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --jobFlavour tomorrow --HCalOnly
python python/submitJobs.py --energy ${energy} --theta 68 --particle kaon0L --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --jobFlavour tomorrow --HCalOnly
# ECal+HCal Endcap
python python/submitJobs.py --energy ${energy} --theta 38 --particle pi- --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --jobFlavour tomorrow
python python/submitJobs.py --energy ${energy} --theta 38 --particle kaon0L --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --jobFlavour tomorrow
# HCal Endcap
python python/submitJobs.py --energy ${energy} --theta 38 --particle pi- --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --jobFlavour tomorrow --HCalOnly
python python/submitJobs.py --energy ${energy} --theta 38 --particle kaon0L --numberOfEvents 10000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --jobFlavour tomorrow --HCalOnly
# ECal+HCal @ 100 GeV
python python/submitJobs.py --energy 100 --theta ${theta} --particle pi- --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim
python python/submitJobs.py --energy 100 --theta ${theta} --particle kaon0L --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim
# HCal @ 100 GeV
python python/submitJobs.py --energy 100 --theta ${theta} --particle e- --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --HCalOnly
python python/submitJobs.py --energy 100 --theta ${theta} --particle pi- --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --HCalOnly
python python/submitJobs.py --energy 100 --theta ${theta} --particle kaon0L --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim --HCalOnly
## Reconstruction
elif [ "$1" == "k4run" ]; then
inputDir="root://128.141.173.81:1094//home/data/FCC/Sim/Resolution"
outputDest="root://128.141.173.81:1094//home/data/FCC/Reco/"
outputDir="Resolution"
energy="2 3 4 5 7 10 13 20 30 50 80 100 150 180"
theta="88 83 78 73 68 63 58 53 48 43 38 33 28 23 18"

# ECal+HCal Barrel | 180 GeV pion is missing
python python/submitJobs.py --energy ${energy} --particle pi- --theta 68 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run
python python/submitJobs.py --energy ${energy} --particle kaon0L --theta 68 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run
# ECal+HCal Endcap
python python/submitJobs.py --energy ${energy} --particle pi- --theta 38 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run
python python/submitJobs.py --energy ${energy} --particle kaon0L --theta 38 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run
# HCal Barrel
python python/submitJobs.py --energy ${energy} --particle pi- --theta 68 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
python python/submitJobs.py --energy ${energy} --particle kaon0L --theta 68 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
python python/submitJobs.py --energy ${energy} --particle e- --theta 68 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
# HCal Endcap
python python/submitJobs.py --energy ${energy} --particle pi- --theta 38 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
python python/submitJobs.py --energy ${energy} --particle kaon0L --theta 38 --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
# ECal+HCal @ 100 GeV
python python/submitJobs.py --energy 100 --theta ${theta} --particle pi- --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run
python python/submitJobs.py --energy 100 --theta ${theta} --particle kaon0L --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run
# HCal @ 100 GeV
python python/submitJobs.py --energy 100 --theta ${theta} --particle pi- --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
python python/submitJobs.py --energy 100 --theta ${theta} --particle kaon0L --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
python python/submitJobs.py --energy 100 --theta ${theta} --particle e- --inputDir ${inputDir} --outputDir ${outputDir} --outputDest ${outputDest} --run k4run --HCalOnly
else
echo "please chose ddsim or k4run"
fi
