#!/usr/bin/env python

# Neighborhood Correlation installation script
# (see http://www.neighborhoodcorrelation.org)

# (C) 2009 Jacob Joseph <jmjoseph@andrew.cmu.edu>,
#          and Carnegie Mellon University
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


from distutils.core import setup, Extension

import sys, os

def get_numpy_include():

    try:
        import numpy
    except ImportError:
        m =  """ERROR: Numpy not found.  This package is required, and may be
downloaded from http://numpy.scipy.org/."""
        print(m, file=sys.stderr)

        return None

    include_dir = numpy.get_include()

    if not os.path.exists( os.path.join(include_dir, 'numpy/arrayobject.h')):
        m = """Numpy appears to be installed, but the reported include path,
'%s' does not contain the header 'numpy/arrayobject.h'""" % include_dir
        print(m, file=sys.stderr)
        return None

    return include_dir


if __name__ == "__main__":

    numpy_include = get_numpy_include()
    
    if numpy_include is None: sys.exit(1)

    setup(name="Neighborhood Correlation",
          version="2.0",
          url="http://www.neighborhoodcorrelation.org",
          author="Jacob Joseph",
          author_email="jmjoseph@andrew.cmu.edu",
          license="GNU General Public License, Version 3 or later",
          packages=['nc'],
          ext_modules=[Extension("nc.nchelparr", ['src/nchelparr.c'])],
          include_dirs = [numpy_include],
          scripts=['NC_standalone']
          )
