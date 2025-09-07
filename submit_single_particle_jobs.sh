## Simulation
if [ "$1" == "ddsim" ]; then
outputDest="root://128.141.173.81:1094//home/data/FCC/Sim/"
outputDir="SingleParticle"
energy="5 10 20 30 50"

python python/submitJobs.py --energy ${energy} --theta 68 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --xrd
python python/submitJobs.py --energy ${energy} --theta 68 --particle kaon0L --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --xrd
python python/submitJobs.py --energy ${energy} --theta 58 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --xrd
python python/submitJobs.py --energy ${energy} --theta 78 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --xrd
python python/submitJobs.py --energy ${energy} --theta 88 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --xrd

## Reconstruction
elif [ "$1" == "k4run" ]; then
inputDir="root://128.141.173.81:1094//home/data/FCC/Sim/SingleParticle/"
outputDest="root://128.141.173.81:1094//home/data/FCC/Reco/"
outputDir="SingleParticle"
energy="5 10 20 30 50"

python python/submitJobs.py --inputDir ${inputDir} --energy ${energy} --theta 68 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --xrd --pandora
python python/submitJobs.py --inputDir ${inputDir} --energy ${energy} --theta 68 --particle kaon0L --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --xrd --pandora
python python/submitJobs.py --inputDir ${inputDir} --energy ${energy} --theta 58 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --xrd --pandora
python python/submitJobs.py --inputDir ${inputDir} --energy ${energy} --theta 78 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --xrd --pandora
python python/submitJobs.py --inputDir ${inputDir} --energy ${energy} --theta 88 --particle pi- --numberOfEvents 5000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --xrd --pandora
else
echo "please chose ddsim or k4run"
fi
