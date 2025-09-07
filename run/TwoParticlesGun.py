#!/usr/bin/python
import argparse
import podio
import podio.root_io
import edm4hep
import math
import numpy as np
from particle import Particle
import os


PDG = {'pi+':211, 'pi-':-211, 'kaon0L':130, 'kaon0S':310, 'kaon0':311}


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--output", dest='output', default="TwoParticlesEvent.root", required=False, help="Output file name")
parser.add_argument("--particles", dest='particles', nargs='+', type=str, default=['kaon0L','pi-'], required=False, help="Particles (pi-, kaon0L)")
parser.add_argument("--energies", dest='energies', nargs='+', type=float, default=[10, 10], required=False, help='Energies of particles in GeV')
parser.add_argument("--thetaMin", dest='thetaMin', type=float, default=0., required=False, help='Min theta angle of the first particle in deg')
parser.add_argument("--thetaMax", dest='thetaMax', type=float, default=0., required=False, help='Max theta angle of the first particle in deg')
parser.add_argument("--alpha", dest='opening_angle', type=float, default=10., required=False, help='Opening angle between the two particles in deg')
parser.add_argument("--nevents", dest='nevents', type=int, default=100, required=True, help="Number of events")
parser.add_argument("--seed", dest='seed', type=int, default=1234, required=False, help="Random seed")

args = parser.parse_args()
np.random.seed(args.seed)

# exit if 2 particles are not provided
if len(args.particles)!=2:
  print('Number of particles is not 2!')
  sys.exit()


def spherical_to_cartesian(phi, theta):
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

def cartesian_to_spherical(vec):
    x, y, z = vec
    theta = np.arccos(z)
    phi = np.arctan2(y, x)
    return phi, theta

def generate_second_particle(phi1, theta1, opening_angle):
    # Convert first particle to Cartesian coordinates
    v1 = spherical_to_cartesian(phi1, theta1)

    # generate a random vector not collinear with v1
    a = spherical_to_cartesian(np.random.uniform(0,2*np.pi), np.random.uniform(0,2*np.pi))

    # check if v1 is close to a, choose another vector to avoid zero vector after projection
    while np.allclose(v1, a):
       a = spherical_to_cartesian(np.random.uniform(0,2*np.pi), np.random.uniform(0,2*np.pi))

    # compute orthogonal vector u to v1
    u = a - np.dot(a, v1) * v1
    # scale it to make unit vector
    u /= np.linalg.norm(u)

    # compute second vector v2 with the given opening angle
    v2 = np.cos(opening_angle) * v1 + np.sin(opening_angle) * u

    # convert back to spherical coordinates
    phi2, theta2 = cartesian_to_spherical(v2)

    return phi2, theta2



# --- Configuration ---
mc_particle_collection_name = "MCParticles"

print("--- Correcting podio file for edm4hep Collections ---")

# Create a podio writer for the output file
writer = podio.root_io.Writer(args.output)

for event_counter in range(args.nevents):
    if event_counter%500==0:
      print("Generating event ",event_counter)

    collection = edm4hep.MCParticleCollection()

    # generate first particle
    pdg = PDG[args.particles[0]]
    part1 = collection.create()
    part1.setPDG(pdg)
    part1.setCharge(Particle.from_pdgid(pdg).charge)
    theta1 = args.thetaMin * np.pi/180
    if args.thetaMin!=args.thetaMax:
       theta1 = np.random.uniform(args.thetaMin,args.thetaMax) * np.pi/180
    phi1 = np.random.uniform(0.,2.*np.pi)
    m = Particle.from_pdgid(pdg).mass/1000.
    e = args.energies[0]
    p = math.sqrt(e**2 - m**2)
    pz = p * math.cos(theta1)
    pt = p * math.sin(theta1)
    px = pt * math.cos(phi1)
    py = pt * math.sin(phi1)
    mom = edm4hep.Vector3d(px,py,pz)
    part1.setMomentum(mom)
    part1.setMass(m)
    part1.setGeneratorStatus(1)

    # generate second particle
    phi2, theta2 = generate_second_particle(phi1, theta1, args.opening_angle * np.pi/180)
    pdg = PDG[args.particles[1]]
    part2 = collection.create()
    part2.setPDG(pdg)
    part2.setCharge(Particle.from_pdgid(pdg).charge)
    m = Particle.from_pdgid(pdg).mass/1000.
    e = args.energies[1]
    p = math.sqrt(e**2 - m**2)
    pz = p * math.cos(theta2)
    pt = p * math.sin(theta2)
    px = pt * math.cos(phi2)
    py = pt * math.sin(phi2)
    mom = edm4hep.Vector3d(px,py,pz)
    part2.setMomentum(mom)
    part2.setMass(m)
    part2.setGeneratorStatus(1)

    R = 2810.5 # mm
    v1 = np.array([R*np.cos(phi1),R*np.sin(phi1),R*np.cos(theta1)/np.sin(theta1)])
    v2 = np.array([R*np.cos(phi2),R*np.sin(phi2),R*np.cos(theta2)/np.sin(theta2)])
    #print(f'Distance at R = {R} mm is {np.linalg.norm(v1-v2)} mm')
    if event_counter%500==0:
       print(f'{event_counter} events are generated...')

    # Create a new frame
    event = podio.Frame()
    # Put the collection into the new frame
    event.put(collection, "MCParticles")

    # Write the entire frame to the output file.
    # The writer automatically handles the data layout.
    writer.write_frame(event,"events")
