
# TOPLOT: Program for plotting and analysis of protein topology.
TOPLOT translates protein structures to topology strings.
The procedure is as follows:
- Read the backbone of the input protein structure.
	Backbone breaks are flagged up.
- Calculate the secondary structure based on Phi/Psi angles.
- Define secondary structure segments based on contiguous stretches
	of sheet (minimum length 4) and helix (minimum length 5).
- Derive the angle between the axes of all pairs of segments that are
	proximate in sequence, i.e. between all segment pairs n-1 and n.
	Antiparallel orientation yields angle 0, parallel oritentation
	yields angle 180. Rotation sense is ignored - all angles are
	positive.
- Define the contact between segment pairs. At least three backbone or Calpha
	atoms need to be closer than the following thresholds (in Angstrom):
```
alpha/alpha 10.5
alpha/beta 8.0
beta/beta 5.5
```
A contact is a binary property.

- Translate the state of each segment pair into a topology character.
	The following alphabet is used (in angle degrees):
```
beta/beta  : A [0-60[, B [60-120[, C [120-180[
alpha/beta : D [0-60[, E [60-120[, F [120-180[
alpha/alpha: G [0-60[, H [60-120[, I [120-180[
```
Segment pairs in contact are upper case as above, those without contact 
are tranformed to lower case.

The topology character of the first segment is set equal to the second 
(owing to the non-existent n-1 segment) but always in lower case.

- Output is the topology string in FASTT format, which is the FASTA format
	containing a topology string. The ouput file has the name of input file
	appended by a '.fastt' suffix.


## Install / Uninstall
Please read the general 'INSTALL' instructions.
Documentation is provided in this README file. To generate formatted and
technical documentation, execute 'doxygen doxygen.cfg' in the 'src' directory.
Documentation files for HTML, LaTeX and man pages will be created in the 
'doc' subdirectory. The latex documentation is completed by executing 
'make pdf' in the 'doc/latex' directory ('refman.ps', 'refman.pdf').


## Usage
```
toplot [--pdb ...] [OPTIONS ...]
  OPTIONS
	--help
```

- OUTPUT: Topology string of protein


## Exit Code and Output Streams
- 0 : Clean termination.
- 1 : Any substantial error


## Copyright and Authors
(C) 2006-2017 Jens Kleinjung


## Availability and License
### Availability
The program is made available under the GNU Public License for academic
scientific purposes, under the condition that proper acknowledgement 
is made to the authors of the program in publications resulting from the use 
of the program.

### License
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

