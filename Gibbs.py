#!/usr/bin/env python3
""" Performs multiple rounds of Gibbs Sampling to find best profile
of k-mers in a file containing DNA sequences"""

import random
import mymotif as mymotif

def Gibbs(filename, k, max_iterate):
    """Perform one round of Gibbs Sampling to find consensus sequence of
    k-mer motifs in a file of DNA strins"""

    f = open(filename)
    Dna = [line.rstrip("\n") for line in f]
    #print(Dna)

    # create list of random k-mers for each dna string
    motifs = []
    for i in range(len(Dna)):
        dna = Dna[i]
        j = random.randint(0,len(dna)-k)
        motif = dna[j:j+k]
        motifs.append(motif)
            
    # initializations for the for-loop
    last_profile = mymotif.profile(motifs)
    last_motifs = motifs
    last_H = mymotif.entropy(last_motifs)
    
    iter_num = 0
    while iter_num < max_iterate:
        iter_num += 1
        # select random i from [0:len(Dna)]
        i = random.randint(0, len(Dna)-k)
        # motifs2 is motifs with motifs[i] missing
        motifs2 = motifs[0:i] + motifs[i+1:]
        # create profile (derived from counts matrix + pseudo counts) of motifs2
        profile2 = mymotif.profile(motifs2,1)
        # get substrings of Dna[i]
        substr_list = mymotif.substrings(Dna[i],k)
        
        
        # get list of probabilities for the substrings
        probs = []
        for substr in substr_list:
            prob = mymotif.prob_motif(substr, profile2)
            probs.append(prob)
            sum_probs = sum(probs)
        for i in range(len(probs)):
            probs[i] = probs[i]/sum_probs
        # randomly select one of these substrings
        new_motif = mymotif.select_new(substr_list, probs)
        # replace motifs[i] with new_motif
        motifs = motifs[0:i] + [new_motif] + motifs[i+1:]
        # Calculate informational entropy (H)
        profile1 = mymotif.profile(motifs)
        H = mymotif.entropy(profile1)
        # stop if converged or too many iterations
        if H == last_H or iter_num == max_iterate:
            return (H, motifs, iter_num)
        last_H = H

def Gibbs_N(filename, k, max_iterate, num_trials):
    best_H, best_motifs, best_iter_num = Gibbs(filename, k, max_iterate)
    for trial in range(num_trials):
        H, motifs, iter_num = Gibbs(filename, k, max_iterate)
        if H < best_H:
            best_H = H
            best_motifs = motifs
            best_iter_num = iter_num
    return(best_iter_num, best_H, best_motifs)


if __name__ == "__main__":
    filename = "dna_seqs.txt"
    k= 6
    max_iterate = 100
    num_trials = 100
    best_iter_num, best_H, best_motifs = Gibbs_N(filename, k, max_iterate, num_trials)

    print("Number of trials =", num_trials)
    print("Best one converged at step #", best_iter_num)
    print("Entropy:", round(best_H,3))
    print("Consensus:", mymotif.consensus(best_motifs))
    best_profile = mymotif.profile(best_motifs)
    mymotif.print_matrix(best_profile)
