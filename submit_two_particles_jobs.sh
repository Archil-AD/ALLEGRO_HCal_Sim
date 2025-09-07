## Simulation
if [ "$1" == "ddsim" ]; then
outputDest="root://128.141.173.81:1094//home/data/FCC/Sim/"
outputDir="TwoParticlesEvent"
theta="88 83 78 73 68 63 58 53"
python python/submitJobs.py --alpha 2 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticlesEvent --xrd
python python/submitJobs.py --alpha 5 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticlesEvent --xrd
python python/submitJobs.py --alpha 10 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticlesEvent --xrd
python python/submitJobs.py --alpha 15 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticlesEvent --xrd
python python/submitJobs.py --alpha 25 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticlesEvent --xrd
python python/submitJobs.py --alpha 35 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run ddsim  --jobFlavour tomorrow --twoParticlesEvent --xrd

## Reconstruction
elif [ "$1" == "k4run" ]; then
inputDir="root://128.141.173.81:1094//home/data/FCC/Sim/TwoParticlesEvent/"
outputDest="root://128.141.173.81:1094//home/data/FCC/Reco/"
outputDir="TwoParticlesEvent"
python python/submitJobs.py --inputDir $inputDir --alpha 2 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticlesEvent --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 5 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticlesEvent --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 10 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticlesEvent --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 15 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticlesEvent --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 25 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticlesEvent --xrd --jobFlavour longlunch
python python/submitJobs.py --inputDir $inputDir --alpha 35 --energy 10 30 --theta $theta --particle kaon0L --numberOfEvents 2000 --outputDest ${outputDest} --outputDir ${outputDir} --run k4run  --jobFlavour tomorrow --twoParticlesEvent --xrd --jobFlavour longlunch
else
echo "please chose ddsim or k4run"
fi
