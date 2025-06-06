# Changelog

## [1.0.0] - 2025-06-03
### Added
- Complementary tool "curvaslarge" (sources and compilation). This is a modification of the original Fernando Flores version, in which more memory is allocated in order to work with larger data arrays (higher sampling rates in the output).

## [0.0.0] - 2024-09-13
First version, as described in Mauro S. Maza's doctoral thesis.

### Main features
Intended to simulate 3-bladed LHAWTs, featruring Fernando G. Flores' Simpact FEM software for the structure response, and Cristian Gebhardt's UVLM implementation for aerodynamic load calculation.
Tower is represented with beam elements, while nacelle and hub are represented as concentrated mases (rigid bodies). Blades might be analized with beam elements (GEBT) or with 3-noded shell elements without rotational DoFs.
