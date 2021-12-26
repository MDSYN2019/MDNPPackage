import argparse
from . import __version__
import logging

from . import solver
from .topology import gen_molecule_smi, gen_molecule_sdf

import sys

def checkArgs(args):

  if not args.sdf and not args.smi:
    parser.error("run requires --sdf or --smi")

  if not args.molname:
    parser.error("run requires --mol")

  if not args.topfname:
    parser.error("run requires --top")

parser = argparse.ArgumentParser(prog='auto_martini', description='Generates Martini force field for atomistic structures of small organic molecules',
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='''Developers:\n===========\nTristan Bereau (bereau [at] mpip-mainz.mpg.de)\nKiran Kanekal (kanekal [at] mpip-mainz.mpg.de)
Andrew Abi-Mansour (andrew.gaam [at] gmail.com)''')

