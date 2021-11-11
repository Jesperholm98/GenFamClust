#!/usr/bin/env python

# Neighborhood Correlation common routines
# (see http://www.neighborhoodcorrelation.org)

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

import sys, time, datetime
from math import log, log10
import numpy

class nc_base( object):
    """Neighborhood Correlation routines"""

    blast_info = None          # dict of BLAST run parameters
    nc_info = None             # dict of NC run parameters    
    
    var_cache = None           # Sequence variance: {'seq_id': var_x, ...}
    E_cache = None             # {'seq_id': (n_x, sum_x, E_x), ...}

    def __init__(self):

        self.blast_info = {}
        self.nc_info = {}

        self.var_cache = {}
        self.E_cache = {}

        # try to compile and use the C module.  Otherwise, fail back
        # to python methods
        try:
            from . import nchelparr
            self.calculate_nc_inner = nchelparr.calculate_nc_inner
        except:
            m = """WARNING(IN NC_BASE): Neighborhood Correlation C helper not found.  Using the pure
Python implementation will be significantly (~90 times) slower.  See
the README documentation for instructions on building."""
            print(m, file=sys.stderr)
            from . import nchelparr_python
            self.calculate_nc_inner = nchelparr_python.calculate_nc_inner
        
            
    def new_run( self, br_id, e_thresh, nc_thresh,
                 smin_factor=0.95, self_hits=1, num_residues=None,
                 score_type='bit_score'):

        self.nc_info = { 'br_id' : br_id,
                         'e_thresh' : e_thresh,
                         'nc_thresh' : nc_thresh,
                         'smin_factor' : smin_factor,
                         'score_type' : score_type,
                         'self_hits' : self_hits}

        if num_residues is not None:
            self.blast_info['num_residues'] = num_residues

        
    def calculate_smin( self):

        assert len(self.blast_info)>0, """blast_info not initialized.
fetch_score_dict_flatfile() must first be executed, or a database call made."""

        # Fetch blast run information
        num_residues = self.blast_info['num_residues']
        blast_e_param = self.blast_info['params']['expectation']
        smin_factor = self.nc_info['smin_factor']

        smin = (log( (num_residues ** 2) / blast_e_param)
                / log(2))

        # reduce the true smin by some factor
        smin = smin * smin_factor
        self.nc_info['smin'] = smin
        self.nc_info['logsmin'] = log10( smin)


    def E_x(self, score_dict, seq_id):
        "Calculate and cache n_x, sum_x, and E_x"

        if seq_id in self.E_cache:
            return self.E_cache[seq_id]

        if not 'smin' in self.nc_info: self.calculate_smin()
        scores = score_dict[seq_id][1]

        num_seqs = self.blast_info['num_sequences']
        logsmin = self.nc_info['logsmin']

        n_x = scores.size
        others = num_seqs - n_x
        sum_x = numpy.sum(scores)
        E_x = (sum_x + logsmin * others) / num_seqs

        self.E_cache[ seq_id] = (n_x, sum_x, E_x)
        return (n_x, sum_x, E_x)


    def var_x(self, score_dict, seq_id):
        "Calculate and cache the variance of one sequence."
        
        if seq_id in self.var_cache:
            return self.var_cache[seq_id]

        # Fetch the mean
        (n_x, sum_x, E_x) = self.E_x(score_dict, seq_id)

        scores = score_dict[seq_id][1]
        num_seqs = self.blast_info['num_sequences']
        logsmin = self.nc_info['logsmin']

        n_others = num_seqs - n_x

        sum_squared = numpy.sum( (scores - E_x)**2 )
        sum_squared += ((logsmin - E_x)**2) * n_others
        var = sum_squared / (num_seqs - 1)

        # Correct for precision issues, or if no neighbors exist
        if var < 0 or n_x==0:
            var = 0

        self.var_cache[seq_id] = var
        return var

    def cache_e_var(self, score_dict, query_seqs=None):
        self.E_cache = {}
        self.var_cache = {}

        if query_seqs is None: query_seqs = list(score_dict.keys())
        for seq_id in query_seqs:
            self.E_x( score_dict, seq_id)
            self.var_x( score_dict, seq_id)
        return (self.E_cache, self.var_cache)


    def calculate_nc(self, score_dict, query_seqs=None,
                     target_seqs=None, callback=None,
                     fd_out=None, seq_id_map=None):
        """Calculate NC scores for all pairs between query_seqs and
target_seqs, using the BLAST scores in score_dict.

* score_dict should be of the form:
  {'seq_id': (numpy.array([id0, id1, ...], dtype=int32),
              numpy.array([score0, score1, ...], dtype=float64)),
   ... }

  All seq_ids are expected to be of type int32.  The arrays should be
  sorted in order of ascending id.

Optional arguments:

* query_seqs and target_seqs should be lists of int32 seq_ids.  When
  unspecified, query_seqs or target_seqs will default to all keys in
  score_dict.

* callback(seq_id_0, n_0, [(seq_id_1, n_1, n_01, nc_score), ...]) will
  be called for each query sequence.  n_0 and n_1 are the number of
  BLAST hits unique to seq_id_0 and seq_id_1, respectively.  n_01 is
  the size of the common neighborhood.

* fd_out will be written with 'seq_id_0 seq_id_1 score\n'.  This must
  be an open file object.  This function will not close the file.

* seq_id_map should be a dictionary mapping from (int32) seq_id to
  some other identifier, of any type with a __str__() method.  This
  mapping will be applied before writing to fd_out."""

        if query_seqs is None:
            query_seqs = list(score_dict.keys())

        # all sequences that will be accessed should have data cached.
        # If both the queries and targets are not specified, cache all
        # sequences.
        if query_seqs is not None and target_seqs is not None:
            cache_seqs = set(query_seqs).union( target_seqs)
        else:
            cache_seqs = list(score_dict.keys())

        #print "query_seqs:", query_seqs
        #print "target_seqs:", target_seqs
        #print "cache_seqs:", cache_seqs

        if target_seqs is not None:
            target_seqs = set( target_seqs)

        # cache all sequence expectations, variances:
        (E_cache, var_cache) = self.cache_e_var(
            score_dict,
            query_seqs = cache_seqs)

        print("       # Queries:", len(query_seqs), file=sys.stderr)
        print("       # Targets:", len(target_seqs) if target_seqs is not None else len(score_dict), file=sys.stderr)
        print("    NC Threshold:", self.nc_info['nc_thresh'], file=sys.stderr)
        print("            Smin:", self.nc_info['smin'], file=sys.stderr)
        print("Sequence DB size:", self.blast_info['num_sequences'], file=sys.stderr)

        t = time.time()
        for i,seq_id_0 in enumerate(query_seqs):
            if (i > 0 and i < 1000 and i % 100==0) or (i > 0 and i % 1000==0):
                print("Queries performed: %d, Elapsed time: %s" % (
                    i, datetime.timedelta(seconds = time.time() - t)))

            if not seq_id_0 in score_dict:
                print("Query '%d' not found in score_dict" % seq_id_0, file=sys.stderr)
                continue
            
            (n_0, nc_hitlist) = self.calculate_nc_inner(
                score_dict, E_cache, var_cache,
                self.nc_info, self.blast_info,
                seq_id_0, 
                target_seqs)

            if fd_out is not None:
                for (seq_id_1, n_1, n_01, nc_score) in nc_hitlist:
                    print("%s %s %r" % (
                        seq_id_map[seq_id_0] if seq_id_map is not None else seq_id_0,
                        seq_id_map[seq_id_1] if seq_id_map is not None else seq_id_1,
                        nc_score), file=fd_out)

            if callback is not None:
                callback( seq_id_0, n_0, nc_hitlist)

        print("Queries performed: %d, Elapsed time: %s" % (
            i, datetime.timedelta(seconds = time.time() - t)))

