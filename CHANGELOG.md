# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.2.0] - 2024-05-14

MINOR release of SURFMAP - minor breaking changes

### Added
- command line option to set min/max values for the bfactor color scale (`--bfactor-min-value`, `--bfactor-max-value`)

### Changed
- change command line option to set max absolute value for electrostatics color scale [BREAKING CHANGE]
  - from `--color-max-val` to  `--elec-max-value`


## [2.1.0] - 2023-12-27

MINOR release of SURFMAP - no breaking changes

### Added
- command line option to set min/max values for the electrostatics color scale

### Patched
- change the Proline value in the Kyte and Doolittle scale (from 1.6 to -1.6)


## [2.0.0] - 2023-02-27

MAJOR release where SURFMAP has been completely restructured as a python package.

### Added
-  a new docker CLI handler
-  surfmap options `-pqr` and `-ff`
- `surfmap` executable command
- `extract_interface` utility script
- `write_pdb_bs` utility script
- MANIFEST.in
- examples for new utility scripts
- setup.py
- Dockerfile
- check for surfmap requirements (apbs, r)
- optional environment variable for surfmap image version
- logging with 3 verbose levels
- surfmap.log file as output
- multival_csv2pdb as binary

### Changed
- README
- LICENCE from GPL to MIT
- install process
- docker repository name
- pdb_to_xyzr from awk to python3
- copyright notices
- path-related instructions in R scripts
- `compute_shell.sh` from bash to a python function
- electrostatics computation as a python function
- version format

### Removed
- use of `surfmap.py`
- use of `run_surfmap.py`
- electrostatics computation from shell computation
- awk/nawk dependency
- unused JS scripts/packages

### Deprecated
- everything included in lower versions is now deprecated


## [1.5] - 2023-02-27








