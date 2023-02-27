from pathlib import Path
import subprocess
import sys
from typing import List

from surfmap import SURFMAP_DOCKER_VERSION


class DockerCLI:
    DOCKER_REPOSITORY = 'lopesi2bc/surfmap'
    TAG = SURFMAP_DOCKER_VERSION

    def __init__(self, args, input_args: List[str], output_args: str, out_dirname: str, input_dir: str='/home/surfmap/input', output_dir: str='/home/surfmap/output') -> None:
        self.cli = ['docker', 'run', '-it', '--rm']

        self.prog = Path(sys.argv[0]).name
        self.input_args = input_args
        self.output_args = output_args
        self.input_dir = input_dir
        self.output_dir = output_dir

        self.input_files = [ args.__getattribute__(x.replace("-", "")) for x in self.input_args]
        self.output = out_dirname
        Path(self.output).mkdir(exist_ok=True, parents=True)

        self._set_cli()

    def _set_cli(self):
        self.cli += self._mounting_points() + [f"{self.DOCKER_REPOSITORY}:{self.TAG}"] + [self.prog] + self._arguments()

    def _mounting_points(self):
        mnt_point = []
        for input_file in self.input_files:
            mnt_point += ["-v", f"{Path(input_file).resolve()}:{self.input_dir}/{Path(input_file).name}"]
        
        mnt_point += ["-v", f"{Path(self.output).resolve()}:{self.output_dir}:rw"]

        return mnt_point

    def _arguments(self):
        cmd_arguments = []

        # add all surfmap command line arguments except --docker
        for argv in sys.argv[1:]:
            if argv == "--docker":
                continue
            cmd_arguments.append(argv)

        # edit arguments relative to input(s)
        for i, input_arg in enumerate(self.input_args, start=0):
            cmd_arguments[cmd_arguments.index(input_arg) + 1] = f"{self.input_dir}/{Path(self.input_files[i]).name}"

        # edit arguments relative to output
        if self.output_args in cmd_arguments:
            cmd_arguments[cmd_arguments.index(self.output_args) + 1] = self.output_dir
        else:
            cmd_arguments += [self.output_args, self.output_dir]

        return cmd_arguments


    def cmd(self):
        return " ".join(self.cli)

    def show(self):
        print(f"\n{self.cmd()}\n")

    def run(self):
        subprocess.call(self.cli)
