#    sgao: Simple Genetic Algorithms for Optimization
#    R library that contains functions for solving non-linear optimization
#    using (simple) genetic algorithms and some appropriate plots
#    Copyright (C) 2003 Kurt Sys
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


seao.terminology<-function() {
  cat ("population ('pop'): A number of experimental set-ups. This can be one batch \n")
  cat ("   or a list of several batches.\n")
  cat ("generation ('gen'): A number of experimental set-ups run at the same time;\n")
  cat ("   a batch of experiments.")
  cat ("individual ('ind'): One experimental set-up with its own set of  parameter\n")
  cat ("   set\n")
  cat ("genome ('genome'): Set of parameter values or parameter set")
  cat ("gene ('gene'): Parameter that can be changed for each experimental setup.\n")
  cat ("allele ('allele'): Value of a parameter in an experimental setup. Combining\n")
  cat ("   the appropriate parameter values results in an individual with highest\n")
  cat ("   fitness.\n")
  cat ("fitness ('fit'): The output of each experimental setup, i.e. the property \n")
  cat ("   you are trying to maximise by searching an optimal parameter set.\n")
}
