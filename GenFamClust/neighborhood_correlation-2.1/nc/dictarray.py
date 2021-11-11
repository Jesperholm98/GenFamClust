#!/usr/bin/env python

# (C) 2011 Jacob Joseph <jmjoseph@andrew.cmu.edu>
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

import sys, csv
from numpy import *

class sortedarray:
    """Maintain a numpy array of elements, and another of values.
    These are ordered by element, and support fast insertion and
    search.  Insertion of a duplicate element should update the value
    of the element with the max of that existing and inserted."""

    def __init__(self, init_size=10, growth_factor=1.5):
        self.growth_factor = growth_factor
        self.nelements = 0
        self.elem = empty( init_size, dtype=int32)
        self.vals = empty( init_size, dtype=float64)

        # Repetitive imports shouldn't have an impact on performance.
        try:
            from . import nchelparr
            self.resolve_dups = nchelparr.resolvesymmdups
        except:
            from . import nchelparr_python
            self.resolve_dups = nchelparr_python.resolvesymmdups

    def __repr__(self):
        return str((self.elem[:self.nelements],
                self.vals[:self.nelements]))

    def grow(self):
        if self.nelements > 1:
            size = int(ceil(self.nelements * self.growth_factor))
        else:
            size = 1
        self.elem.resize( size, refcheck=False)
        self.vals.resize( size, refcheck=False)
        return

    def to_arraytuple(self):
        """Return the underlying data.  Typically called after
        building the structure, when dynamic resizing is no longer
        necessary."""
        return (self.elem, self.vals)

    def insert_unsorted(self, element, value):
        """Append to the end of the array, in preparation for sorting
        in one pass"""

        if self.elem.size < self.nelements + 1:
            self.grow()

        self.elem[ self.nelements] = element
        self.vals[ self.nelements] = value
        
        self.nelements += 1
        return

    def trim(self):
        """Trim the array memory to exactly fit the number of elements
        present."""
        
        self.elem = resize(self.elem, self.nelements)
        self.vals = resize(self.vals, self.nelements)
        return

    def sort(self):
        """Sort, and remove duplicates."""

        order = argsort( self.elem[:self.nelements])
        self.elem[:self.nelements] = self.elem[ order]
        self.vals[:self.nelements] = self.vals[ order]
        return

    def finish(self):
        """Sort, remove duplicates, trim the array memory to its real
        size"""

        self.sort()
        self.nelements = self.resolve_dups(self.elem, self.vals, self.nelements)
        self.trim()
        return


def fetch_score_dict_flatfile( fname):
    """Build a dictionary/array data structure of symmetric scores
from a file of BLAST scores, in the format of
    'seq_id_0 seq_id_1 bit_score',
one per line.

Output (score_dict, seq_id_map), where score_dict is a dictionary of
score arrays:

    {int0: (numpy.array([int1, int2, ...]),
            numpy.array([score1,score2, ...])),
     ... }

And seq_id_map is dictionary mapping from assigned integer to the
original string in the input file:
    {int0: seq_id_0,
     int1: seq_id_1,
     ... }.

Scores will be transformed by application of a log_10."""
    
    # use the sortedarray class for the hits while reading, then
    # convert the result back to standard numpy arrays

    print("Reading BLAST bit-score file: %s" % fname, file=sys.stderr)
    #print >> sys.stderr, time.ctime()

    try:
        from nc import nchelparr

    except:
        m = """WARNING: Neighborhood Correlation C helper not found.  Using the pure
        Python implementation will be significantly (~90 times) slower.  See
        the readme documentation for instructions on building."""
        print(m, file=sys.stderr)
            
        from nc import nchelparr_python

    fd = open(fname)
    sd = {}

    # sequence id map.
    seq_id_map_index = {}  # integers to keys
    next_seq_id = 0
    seq_id_map = {}        # keys to integers

    # sequence id map.
    seq_id_map_index = {}  # integers to keys
    next_seq_id = 0
    seq_id_map = {}        # keys to integers

    for larr in csv.reader(fd, delimiter=" "):
        try:
            (seq_id_0_txt, seq_id_1_txt, score) = larr
        except:
            print("Flatfile line unparsable: '%s'" % l, file=sys.stderr)
            raise

        try:
            seq_id_0 = seq_id_map_index[seq_id_0_txt]
        except KeyError:
            seq_id_0 = next_seq_id
            seq_id_map_index[seq_id_0_txt] = seq_id_0
            seq_id_map[seq_id_0] = seq_id_0_txt
            sd[seq_id_0] = sortedarray()
                
            next_seq_id += 1

        try:
            seq_id_1 = seq_id_map_index[seq_id_1_txt]
        except KeyError:
            seq_id_1 = next_seq_id
            seq_id_map_index[seq_id_1_txt] = seq_id_1
            seq_id_map[seq_id_1] = seq_id_1_txt
            sd[seq_id_1] = sortedarray()
               
            next_seq_id += 1

        score = log10( float(score))

        sd[seq_id_0].insert_unsorted( seq_id_1, score)
        sd[seq_id_1].insert_unsorted( seq_id_0, score)

    fd.close()

    # convert back to standard numpy arrrays
    for k, sa in sd.items():
        sa.finish()
        sd[k] = sa.to_arraytuple()

    print("Completed symmetric score construction.", file=sys.stderr)
    #print >> sys.stderr, time.ctime()
    return (sd, seq_id_map)

