# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2023-02-27

MAJOR release where SURFMAP has been completely restructured as a python package.

### Added
- added a new docker CLI handler
- added surfmap options `-pqr` and `-ff`
- added `surfmap` executable command
- added `extract_interface` utility script
- added `write_pdb_bs` utility script
- added MANIFEST.in
- added examples for new utility scripts
- added setup.py
- added Dockerfile
- added check for surfmap requirements (apbs, r)
- added optional environment variable for surfmap image version
- added logging with 3 verbose levels
- added surfmap.log file as output
- added multival_csv2pdb as binary

### Changed
- changed README
- changed LICENCE from GPL to MIT
- changed install process
- changed docker repository name
- changed pdb_to_xyzr from awk to python3
- changed copyright notices
- changed path-related instructions in R scripts
- changed `compute_shell.sh` from bash to a python function
- changed electrostatics computation as a python function
- changed version format

### Removed
- removed use of `surfmap.py`
- removed use of `run_surfmap.py`
- removed electrostatics computation from shell computation
- removed awk/nawk dependency
- removed unused JS scripts/packages

### Deprecated
- everything included in lower versions is now deprecated


## [1.5] - 2023-02-27








