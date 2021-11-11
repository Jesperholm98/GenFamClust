#!/usr/bin/env python

# Jacob Joseph

# Python implementation of NC calculation, mirroring those in
# nchelparr.c.  These will be much slower than the C helper, but don't
# require compilation.


#-----COMMENTED AWAY FOR PYTHON 3.X VERSIONS-----
#from IPython.Shell import IPShellEmbed


def resolvesymmdups( elem, vals, nelements):
    i=0
    shift=0

    while i + shift + 1 < nelements:
        #print('loop', i, shift, nelements)
        if elem[i] == elem[i+shift+1]:
            
            if vals[i+shift+1] > vals[i]:
                vals[i] = vals[i+shift+1]
            shift += 1
        else:
            if shift > 0:
                elem[i+1] = elem[i+shift+1]
                vals[i+1] = vals[i+shift+1]
            i += 1
    #print(i, shift, nelements)
    return nelements - shift


def covariance_xy( hits_0, scores_0, expect_tup_0,
                   hits_1, scores_1, expect_tup_1,
                   num_seqs, logsmin):

    (n_0, sum_0, Ex) = expect_tup_0
    (n_1, sum_1, Ey) = expect_tup_1

    n_01 = 0
    s01_sum = 0
    s01_prodsum = 0

    ind_0 = 0
    ind_1 = 0
    len_0 = hits_0.size
    len_1 = hits_1.size
    
    # hit lists are in sorted order
    while ind_0 < len_0 and ind_1 < len_1:
        h0 = hits_0[ind_0]
        h1 = hits_1[ind_1]
        if h0 < h1:
            ind_0 += 1
            continue
        elif h0 > h1:
            ind_1 += 1
            continue

        # elif h0 == h1
        n_01 += 1
        s01_sum += scores_0[ind_0]
        s01_sum += scores_1[ind_1]
        s01_prodsum += scores_0[ind_0] * scores_1[ind_1]

        ind_0 += 1
        ind_1 += 1

    n_others = num_seqs - (n_0 + n_1 - n_01)

    Exy = (s01_prodsum + (logsmin**2) * n_others +
           (sum_0 + sum_1 - s01_sum) * logsmin
           ) / num_seqs

    var_xy =  Exy - Ex * Ey

    return (n_01, var_xy)


def calculate_nc_inner(score_dict, E_cache, var_cache,
                       nc_info, blast_info, seq_id_0, target_seqs):

    logsmin = nc_info['logsmin']
    nc_thresh = nc_info['nc_thresh']
    num_seqs = blast_info['num_sequences']

    (hits_0, scores_0) = score_dict[ seq_id_0]
    var_0 = var_cache[ seq_id_0]
    expect_tup_0 = E_cache[ seq_id_0]

    neighs = set( hits_0)
    nneighs = set()

    retlist = []
    accumulate_nneighs = True
    i = 0
    while True:
        if accumulate_nneighs and i < len(hits_0):
            seq_id_1 = hits_0[i]
            i += 1
        else:
            accumulate_nneighs = False

            # if len(nneighs) < 1: break
            try:
                seq_id_1 = nneighs.pop()
            except KeyError:
                break

        if accumulate_nneighs:
            hits_1 = set(score_dict[ seq_id_1][0])

            #print "hits_1 = ", hits_1
            #print "target_seqs = ", target_seqs
            #print "hits_1.intersection(target_seqs):", hits_1.intersection(target_seqs)

            #ipsh = IPShellEmbed(argv=[])
            #ipsh("In nc calculation:")

            # intersect with the target list
            if target_seqs is not None:
                hits_1.intersection_update(target_seqs)

            #print "hits_1", hits_1
            
            # remove immediate neighbors
            hits_1.difference_update( neighs)
            
            nneighs.update( hits_1)

        # check the target list.  neighs may have had items not in it
        if target_seqs is not None and seq_id_1 not in target_seqs:
            continue

        if not seq_id_1 in score_dict:
            print("Skipping %d - %d. seq_id_1 not found in score dictionary.  Are scores symmetric?" % (
                seq_id_0, seq_id_1))
            continue

        (hits_1, scores_1) = score_dict[ seq_id_1]
        var_1 = var_cache[ seq_id_1]
        expect_tup_1 = E_cache[ seq_id_1]

        (n_01, cov_01) = covariance_xy( hits_0, scores_0, expect_tup_0,
                                        hits_1, scores_1, expect_tup_1,
                                        num_seqs, logsmin)

        if (var_0 == 0 or var_1 == 0):
            if n_01 > 0: nc = 1.0
        else:
            nc = cov_01 / ( var_0 * var_1)**0.5

            if nc > nc_thresh:
                retlist.append( (seq_id_1, expect_tup_1[0], n_01, nc))
        
    return (len(hits_0), retlist)


