/* 
   NC calculation using numpy arrays

   (C) 2011 Jacob Joseph <jmjoseph@andrew.cmu.edu>
            and Carnegie Mellon University

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.
*/

#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numpy/arrayobject.h>

/* Remove duplicate elements from lists of elements and values.  Select
   the greatest value of any duplicate. This is used in symmetric
   score calculation.  See dictarray.py */
PyObject *resolvesymmdups( PyObject *s, PyObject *args) {
  PyObject *pyelem, *pyvals, *pynelements;
  int i = 0, shift = 0, nelements;
  int *elem;
  double *vals;
  
  if (!PyArg_ParseTuple( args, "OOO", &pyelem, &pyvals, &pynelements))
    return NULL;

  elem = (int *)PyArray_DATA( pyelem);
  vals = (double *)PyArray_DATA( pyvals);
  nelements = PyLong_AsLong( pynelements);

  while (i + shift + 1 < nelements) {
    if ( elem[i] == elem[i+shift+1]) {

      if ( isgreater( vals[i+shift+1], vals[i]))
	vals[i] = vals[i+shift+1];
    
      shift += 1;
      continue;
    }
    else {
      if (shift > 0) {
	elem[i+1] = elem[i+shift+1];
	vals[i+1] = vals[i+shift+1];
      }
      i += 1;
    }
  }

  return Py_BuildValue("i", nelements - shift);
}

/* Calculate the covariance between two sequences.  Return values in
   n_01, var_01  */
int covariance_xy ( PyObject *hits_0, PyObject *scores_0, PyObject *expect_tup_0,
		    PyObject *hits_1, PyObject *scores_1, PyObject *expect_tup_1,
		    int num_seqs, double logsmin, int *n_01, double *cov) {
  int n_others, ind_0 = 0, ind_1 = 0;
  double s01_sum = 0, s01_prodsum = 0;

  /* Find the intersection of hits from the two sequences.  Keep a
     count of the size of this intersection, n_01, a sum of these
     scores, s01_sum, and a sum of the products of these scores,
     s01_prodsum */
  *n_01 = 0;
  while ( (ind_0 < PyArray_DIM(hits_0,0)) && (ind_1 < PyArray_DIM(hits_1,0)) ) 
    {
      int h0 = ((int *)PyArray_DATA(hits_0))[ind_0];
      int h1 = ((int *)PyArray_DATA(hits_1))[ind_1];

      if (h0 < h1) {
	ind_0++;
	continue;
      }
      else if (h0 > h1) {
	ind_1++;
	continue;
      }

      /* hit found in both sequences */
      *n_01 += 1;
      s01_sum += ((double *)PyArray_DATA(scores_0))[ind_0];
      s01_sum += ((double *)PyArray_DATA(scores_1))[ind_1];
      s01_prodsum += ((double *)PyArray_DATA(scores_0))[ind_0] *
	((double *)PyArray_DATA(scores_1))[ind_1];
      ind_0 += 1;
      ind_1 += 1;
    }

  /*n_others = num_seqs - (n_0 + n_1 - n_01);*/
  n_others = num_seqs - (PyLong_AsLong( PyTuple_GetItem( expect_tup_0, 0) ) +
			 PyLong_AsLong( PyTuple_GetItem( expect_tup_1, 0) ) -
			 *n_01);

  /* Calculate Exy.  First add the pseudo matches of logsmin in both sequences */
  *cov = s01_prodsum + pow(logsmin,2) * n_others;
        
  /* Now add the factor of only one sequence having a match, and the other logsmin */
  *cov += ( PyFloat_AsDouble( PyTuple_GetItem( expect_tup_0, 1)) +
	   PyFloat_AsDouble( PyTuple_GetItem( expect_tup_1, 1)) -
	   s01_sum) * (double)logsmin;

  /* Exy is normalized by the number of sequences */
  *cov /= num_seqs;

  /* covariance = Exy - Ex*Ey */
  *cov -= PyFloat_AsDouble( PyTuple_GetItem( expect_tup_0, 2)) *
    PyFloat_AsDouble( PyTuple_GetItem( expect_tup_1, 2));

  return 0;
}

/* PyLong_AsLong() doesn't work on numpy.int objects, so build a
   python object first.  This is primarily for debug messages with
   printf.  I'm unsure without further checking why PyArray_PyIntAsInt
   shouldn't handle this, but it segfaults.
*/
long cast_long( PyObject *o) {
  PyObject *tmp;
  long l = -1;
  static char *msg = "an integer is required";

  if (PyLong_Check(o)) {
    l = PyLong_AsLong( o);
    goto finish;
  }
  else if (PyLong_Check(o)) {
    l = PyLong_AsLong(o);
    //l = PyInt_AsLong(o);
    goto finish;
  }
  else {
    tmp = PyNumber_Long( o);
    l = PyLong_AsLong( tmp);
    Py_DECREF(tmp);
  }

 finish:
  if ((l) == -1 && PyErr_Occurred()) {
      PyErr_SetString(PyExc_TypeError, msg);
      return -1;
    }
  return l;
}

/* calculate_nc_inner Return a python list of tuples (seq_id_1, n_0,
   n_1, n_01, nc), where nc >= nc_thresh.  keys is an iterator of
   seq_id_1 identifiers; if not None, only scores to these sequences
   will be calculated and returned.
*/
PyObject *calculate_nc_inner( PyObject *s, PyObject *args) {
  PyObject *seq_id_0, *target_set, *retlist;
  PyObject *hits_dict, *e_dict, *var_dict, *nc_info, *blast_info;
  PyObject *tmp;
  PyObject *hits_0, *scores_0, *expect_tup_0;
  double nc_thresh, var_0, logsmin;
  long num_seqs;

  PyObject *neighs, *nneighs;
  int i,j;
  char accumulate_nneighs;

  if (!PyArg_ParseTuple( args, "OOOOOOO", &hits_dict, &e_dict, &var_dict, 
			 &nc_info, &blast_info, &seq_id_0, &target_set)) {
    return NULL;
  }


  /* Return a list of tuples */
  retlist = PyList_New(0);

  /* Fetch parameters for NC calculation */
  logsmin = PyFloat_AsDouble( PyDict_GetItemString(nc_info, "logsmin"));
  nc_thresh = PyFloat_AsDouble( PyDict_GetItemString(nc_info, "nc_thresh"));
  num_seqs = PyFloat_AsDouble( PyDict_GetItemString(blast_info, "num_sequences"));

  /* Fetch the score dictionaries for seq_id_0 */
  tmp = PyDict_GetItem( hits_dict, seq_id_0);
  hits_0 = PyTuple_GetItem( tmp, 0);
  scores_0 = PyTuple_GetItem( tmp, 1);
  expect_tup_0 = PyDict_GetItem( e_dict, seq_id_0);
  var_0 = PyFloat_AsDouble( PyDict_GetItem(var_dict, seq_id_0));

  /* build a python set of neighbors, for fast(?) lookup.  Here, we
     cannot filter by the target sequences since their next-neighbors
     may be of interest but not directly connected. */
  neighs = PySet_New(NULL);
  for (j=0; j < PyArray_DIM(hits_0,0); j++) {
    tmp = PyLong_FromLong( ((int *)PyArray_DATA(hits_0))[j]);
    PySet_Add( neighs, tmp);
    Py_DECREF(tmp);
  }
  
  if (PyErr_Occurred()) {
    fprintf(stderr, "PyErr occurred.  Do something.\n");
  }

  nneighs = PySet_New(NULL);
  i = 0;
  accumulate_nneighs = 1;
  while (1) {
    PyObject *seq_id_1, *hits_1, *scores_1, *expect_tup_1;    
    int n_01;
    double cov_01, var_1, nc;

    /* Iterate through the immediate neighbors.  Accumulate 'next'
       neighbors as we go.  Once the neighbors are finished, iterate
       through the nneighs, and don't accumulate */
    
    if (accumulate_nneighs && i < PyArray_DIM(hits_0,0)) {
      seq_id_1 = PyLong_FromLong( ((int *)PyArray_DATA(hits_0))[i] );
      i += 1;
    }
    else {
      if (PySet_Size( nneighs) < 1) break;

      accumulate_nneighs = 0;

      seq_id_1 = PySet_Pop( nneighs);
      if (seq_id_1 == NULL) break;  /* done */
    }

    if (!PyMapping_HasKey( hits_dict, seq_id_1)) {
      fprintf(stderr, "Skipping %ld - %ld. seq_id_1 not found in score dictionary.  Are scores symmetric?\n",
	     cast_long(seq_id_0), cast_long(seq_id_1));
      Py_DECREF(seq_id_1);
      continue;
    }

    /* Lookup the hit arrays, and cached Ex, varx for seq_id_1 */
    tmp = PyDict_GetItem( hits_dict, seq_id_1);
    hits_1 = PyTuple_GetItem( tmp, 0);
    scores_1 = PyTuple_GetItem( tmp, 1);

    if (accumulate_nneighs) {
      for (j=0; j < PyArray_DIM(hits_1,0); j++) {
	tmp = PyLong_FromLong( ((int *)PyArray_DATA(hits_1))[j]);

	/* check the target list */
	if ( PyAnySet_Check( target_set) && !PySet_Contains( target_set, tmp)) {
	  goto skip_nneigh;
	}

	/* check that this isn't an immediate neighbor */
	if (! PySet_Contains(neighs, tmp)) {
	  PySet_Add( nneighs, tmp);
	}

      skip_nneigh:
	Py_DECREF(tmp);
      }
    }

    /* check the target list.  neighs may have had items not in it */
    if ( PyAnySet_Check( target_set) && !PySet_Contains( target_set, seq_id_1)) {
      goto skip_nc_calculation;
    }

    expect_tup_1 = PyDict_GetItem( e_dict, seq_id_1);
    var_1 = PyFloat_AsDouble( PyDict_GetItem(var_dict, seq_id_1));
    if (covariance_xy( hits_0, scores_0, expect_tup_0,
		       hits_1, scores_1, expect_tup_1,
		       num_seqs, logsmin, &n_01, &cov_01)) {
      fprintf(stderr, "covariance failed: %ld, %ld\n", 
	      cast_long(seq_id_0),
	      cast_long(seq_id_1));
      Py_DECREF(seq_id_1);
      Py_DECREF(retlist);
      Py_DECREF(neighs);
      Py_DECREF(nneighs);
      return NULL;
    }

    /* Calculate NC. A sequence could happen to have zero variance.
       If n_01 is zero, NC==0.  If the two sequences have zero
       variance, but n_01 > 0, NC=1. */
    if ( 0 == var_0 || 0 == var_1) {
      if (n_01 > 0) nc = 1.0;
      else nc = 0;
    }
    else {
      nc = cov_01 / pow( var_0 * var_1, 0.5);
    }

    /* Only keep those pairs above a threshold */
    if ( isgreaterequal(nc, nc_thresh) ) {
      tmp = Py_BuildValue("(OOid)", seq_id_1,
			  PyTuple_GetItem(expect_tup_1, (Py_ssize_t)0),
			  n_01, nc);
      PyList_Append(retlist, tmp);
      Py_DECREF(tmp);
    }
  
  skip_nc_calculation:
    Py_DECREF(seq_id_1);
  }

  if (PyErr_Occurred()) {
    fprintf(stderr, "PyErr occurred.  Do something.\n");
  }

  tmp = Py_BuildValue("(OO)", 
		      PyTuple_GetItem(expect_tup_0, (Py_ssize_t)0),
		      retlist);
  Py_DECREF(retlist);
  Py_DECREF(neighs);
  Py_DECREF(nneighs); 

  return tmp;
}


static PyMethodDef nchelpmethods [] = {
  {"calculate_nc_inner", calculate_nc_inner, METH_VARARGS},
  {"resolvesymmdups", resolvesymmdups, METH_VARARGS},
  { NULL, NULL}
};

static struct PyModuleDef nchelparrmodule = {
    PyModuleDef_HEAD_INIT,
    "nchelparr",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,      /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    nchelpmethods
};

PyMODINIT_FUNC PyInit_nchelparr(void)
{
    return PyModule_Create(&nchelparrmodule);
};
  
