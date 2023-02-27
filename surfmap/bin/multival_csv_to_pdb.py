#!/usr/bin/env python

import argparse
from pathlib import Path
from typing import Union


def run(ifile: Union[str, Path], ofile: Union[str, Path]):
    with open(ofile, "w+") as outfile:
        with open(ifile) as f:
            for i, line in enumerate(f, start=1):
                lsp = line.split(",")
                numbers = [format(float(x), '.3f') for x in lsp]

                nb = '{0: >7}'.format(str(i))
                x = '{0: >8}'.format(str(numbers[0]))
                y = '{0: >8}'.format(str(numbers[1]))
                z = '{0: >8}'.format(str(numbers[2]))

                if len(numbers) == 4:
                    bf = '{0: >13}'.format(str(numbers[3]))
                else:
                    bf = '{0: >13}'.format("0.000")

                stl = "ATOM" + nb + "  O   POS A   1    " + x + y + z + bf+ "\n"
                outfile.write(stl)


def main():
    parser = argparse.ArgumentParser(description="Script to convert multivalue output csv file to pdb file.")
    parser.add_argument("-i", required=True, help="csv file from multivalue script (APBS suite)")
    parser.add_argument("-o", required=True, help="output pdb file")
    args = parser.parse_args()

    ifile = args.i
    ofile = args.o

    run(ifile=ifile, ofile=ofile)

if __name__ == "__main__":
    main()
