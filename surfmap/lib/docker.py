"""
surfmap -pdb 1g3n.pdb -tomap all 

surfmap -pdb 1g3n.pdb -tomap all --docker

if "--docker":
then map required io mount point between local and container file systems 

docker run --rm -it -v Path(params.pdbarg).parent.resolve():/input Path(params.outdir).resolve():/output nchenche/surfmap:2.0 surfmap -pdb /input/params.pdbname -tomap stickiness -d /output/

"""
from pathlib import Path
import sys

from surfmap.lib.parameters import Parameters


DOCKER_REPOSITORY = 'nchenche/surfmap'
IMAGE_TAG = '2.0'
DOCKER_IMAGE = DOCKER_REPOSITORY + ':' + IMAGE_TAG


def cli_adapter(params: Parameters) -> Parameters:
    docker_cli = [
        'docker',
        'run',
        '-it',
        '--rm',
        '--volume',
        f"{Path(params.pdbarg).parent.resolve()}:/input/:rw",
        '--volume',
        f"{Path(params.outdir).resolve()}:/output/:rw",
        DOCKER_IMAGE,
        'surfmap',
    ]

    for argv in sys.argv[1:]:
        if argv == "--docker":
            continue
        docker_cli.append(argv)

    input_arg = [ x for x in docker_cli if x in ["-pdb", "-mat"] ][0]
    docker_cli[docker_cli.index(input_arg) + 1] = f"/input/{params.pdbname}"

    if "-d" in docker_cli:
        docker_cli[docker_cli.index("-d") + 1] = "/output"
    else:
        docker_cli += ["-d", "/output"]

    return docker_cli