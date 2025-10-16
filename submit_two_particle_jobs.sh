## Simulation
if [ "$1" == "ddsim" ]; then
outputDest="root://pclcd27:1094//home/data/FCC/Sim/"
outputDir="combined"
theta="88 78 68 58"
#theta="48"
python python/submitJobs.py --alpha 5 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticleEvents --xrd
python python/submitJobs.py --alpha 10 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticleEvents --xrd
python python/submitJobs.py --alpha 15 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticleEvents --xrd
python python/submitJobs.py --alpha 25 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticleEvents --xrd

## Reconstruction
elif [ "$1" == "k4run" ]; then
inputDir="root://pclcd27:1094//home/data/FCC/Sim/combined/"
outputDest="root://pclcd27:1094//home/data/FCC/Reco/"
outputDir="combined"
theta="88 78 68 58"
#theta="48"
python python/submitJobs.py --inputDir $inputDir --alpha 5 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticleEvents --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 10 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticleEvents --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 15 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticleEvents --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 25 --energy 30 10 --theta $theta --particle kaon0L --numberOfEvents 1000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticleEvents --xrd --jobFlavour longlunch
else
echo "please chose ddsim or k4run"
fi
