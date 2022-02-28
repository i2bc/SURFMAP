#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, subprocess
import docker

DOCKER_REPOSITORY = 'nchenche/surfmap'
IMAGE_TAG = 'latest'
DOCKER_IMAGE = DOCKER_REPOSITORY + ':' + IMAGE_TAG

OPTIONAL_ARGUMENTS = {
    'proj': '-proj',
    'coords': '-coords',
    'res': '-res',
    'rad': '-rad',
    's': '-s', 
    'nosmooth': '--nosmooth',
    'png': '--png',
    'keep': '--keep'
    }
DATA_VOLUME_PATH = "/home/surfmap/data/"  # volume mount point in docker container

update_parser = argparse.ArgumentParser()
update_parser.add_argument("--update", required=False, action="store_true", help = "Update your docker image version")
update_args = update_parser.parse_args()

if update_args.update:
    try:
        client = docker.client.from_env()
    except:
        print('Error, it seems like the docker daemon service is not running...')
        print('If you have docker installed on your machine, you can start your docker service as follows:')
        print(' - on linux system type in  terminal "sudo service docker start"')
        print(' - on a Windows or macOS system, simply open your docker desktop app GUI')

    try:
        print('Updating your local SURFMAP docker image with the latest availabe in {}'.format(DOCKER_IMAGE))
        client.images.pull(repository=DOCKER_REPOSITORY, tag=IMAGE_TAG)
        print('Updating process done. You are now using the latest version of SURFMAP\n')
    except:
        print('Error while trying to pull the latest docker image in {}\n'.format(DOCKER_IMAGE))
    
    exit()

parser = argparse.ArgumentParser()
parser.add_argument("-pdb",required = True, help = "Input pdb file (path + file name)")
parser.add_argument("-tomap", type = str, required = True, choices = set(("all", "stickiness", "kyte_doolittle", "wimley_white", "electrostatics", "circular_variance", "circular_variance_atom", "bfactor", "binding_sites")), help = "Choice of the scale. Argument must be one of the following: stickiness; kyte_doolittle; wimley_white; electrostatics; circular_variance; bfactor; binding_sites; all")
parser.add_argument("-proj", type = str, required = False, choices = set(("sin", "moll", "aitoff", "cyl")), help = "Choice of the projection. Argument must be one of the following: sinusoidal; mollweide; aitoff; cylequalarea")
parser.add_argument("-coords", help = argparse.SUPPRESS)
parser.add_argument("-res", help = "File containing a list of residues to map on the projection. Format must be the following: col 1 = chain id; col 2 = res number; col 3 = res type")
parser.add_argument("-rad", required = False, help = "Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the smoother the surface (default: 3.0)")
parser.add_argument("-d", required = False, help = "Output directory where all files will be written (default: './output_SURFMAP_$pdb_$tomap' where $pdb and $tomap are the inputs given to -pdb and -tomap arguments, respectiveley)")
parser.add_argument("-s", required = False, help = "Size of a grid cell. Necessary that 180%%cellsize == 0 (default: 5.0)")
parser.add_argument("--nosmooth", action = "store_true", help = "If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!)")
parser.add_argument("--png", action = "store_true", help = "If chosen, a map in png format is computed (default: only pdf format is generated)")
parser.add_argument("--keep", action = "store_true", help = "If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept)")
args = parser.parse_args()

pdbarg = args.pdb
tomap = args.tomap

docker_template_cmd = [
    'docker',
    'run',
    '-it',
    '--rm',
    '--volume',
    os.path.dirname(os.path.abspath(args.pdb)) + ':' + '/input/' + ':rw',
    '--volume',
    os.getcwd() + ':' + '/output/' + ':rw',
    DOCKER_IMAGE,
    '-pdb',
    '/input/' + os.path.basename(args.pdb),
    '-tomap',
    args.tomap,
    '-d',
    '/output/' + 'output_SURFMAP_' + os.path.basename(args.pdb).split('.pdb')[:-1][0] + '_' + args.tomap
]

if args.d:
    docker_template_cmd[7] = os.path.dirname(os.path.abspath(args.d)) + ':' + '/output/' + ':rw'
    docker_template_cmd[14] = '/output/' + os.path.basename(args.d)

for opt_args in OPTIONAL_ARGUMENTS:
    if opt_args in args.__dict__ and args.__dict__.get(opt_args):
        docker_template_cmd.append(OPTIONAL_ARGUMENTS[opt_args])
        if opt_args not in ['nosmooth', 'png', 'keep']:
            if opt_args == 'res':
                docker_template_cmd.append('/input/' + args.__dict__.get(opt_args))
            else:
                docker_template_cmd.append(args.__dict__.get(opt_args))



# Run docker container of surfmap
print('\n' + ' '.join(docker_template_cmd) + '\n')
subprocess.call(docker_template_cmd)
