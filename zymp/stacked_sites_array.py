import numpy as np
from dnachisel.biotools import reverse_complement
from Bio import Restriction
import proglog

from .tools import (enzymes_to_dna_pattern, get_enzymes_ATGC_sequences,
                    find_patterns_matches, reverse_complement)

def sites_difference(site1, site2):
    """Return minimal sequence of nucleotides that should be added at the end
    of site1 to make site2 appear.""" 
    for i in range(len(site2), -1, -1):
        if site2[:i] == site1[-i:]:
            return site2[i:]
    return site2

def _enzymes_names_to_distances_graph(enzymes_names):
    enzymes_names = enzymes_names[:]
    np.random.shuffle(enzymes_names)
    enzymes_sites = get_enzymes_ATGC_sequences(enzymes_names)
    patterns = enzymes_to_dna_pattern(enzymes_names)
    core_enzymes = {}
    for e1 in list(enzymes_names):
        site1 = enzymes_sites[e1]
        rev_site1 = reverse_complement(site1)
        for e2 in list(enzymes_names):
            if e1 == e2:
                continue
            site2 = enzymes_sites[e2]
            pattern2 = patterns[e2]
            if any([site1 == site2,
                    site2 in site1,
                    rev_site1 == site2,
                    site2 in rev_site1,
                    len(pattern2.find_matches(site1)),
                    len(pattern2.find_matches(rev_site1)),
                    ]):
                if e1 not in core_enzymes:
                    core_enzymes[e1] = []
                if e2 not in core_enzymes[e1]:
                    core_enzymes[e1].append(e2)
    graph = {}
    for e1 in enzymes_names:
        site1 = enzymes_sites[e1]
        for e2 in enzymes_names:
            if e1 == e2:
                continue
            site2 = enzymes_sites[e2]
            diff = sites_difference(site1, site2)
            graph[(e1, e2)] = dict(diff=diff, dist=len(diff))
    return graph, enzymes_sites, core_enzymes

def _enzyme_path_to_sequence(path, graph, enzymes_sites):
    """Converts a path of successive enzymes into a sequence""" 
    return "".join([enzymes_sites[path[0]]] + [
        graph[(n1, n2)]['diff']
        for n1, n2 in zip(path, path[1:])
    ])

class NoSuitableSequenceFound(Exception):
    pass

def stacked_sites_array(enzymes_names, forbidden_enzymes=(), unique_sites=True,
                        tries=10, logger='bar', success_condition=None):
    """Generate a sequence

    Parameters
    ----------

    enzymes_names
      Names of enzyme sites which the algorithm should try to include in the
      sequence.
    
    forbidden_enzymes
      List of enzymes names whose sites should not be in the sequence.
    
    unique_sites
      If True all the enzymes in enzymes_names will have no more than one site
      in the final sequence.
    
    
    tries
      Number of tries. the sequence returned is the one with the best score,
      i.e. with the least enzymes left over, or the shortest sequence in case
      of equality.
      
    success_condition
      A function evaluated at the end of each try to validate the designed
      sequence. Should be of the form ``f(seq, inseq, notinseq) => True/False``
      where ``seq`` is the final sequence, ``inseq`` the set of enzymes in the
      sequence, and ``notinseq`` the set of enzymes not in the sequence.
    
    logger
      Either "bar" for a progress bar, None for no progress logger, or any
      Proglog progress logger.
    
    Returns
    -------

    A triplet (sequence, enzymes_in_sequence, enzymes_not_in_sequence)
    """

    enzymes_names = sorted(set(enzymes_names))
    logger = proglog.default_bar_logger(logger, min_time_interval=0.2)
    patterns = enzymes_to_dna_pattern(list(enzymes_names) +
                                      list(forbidden_enzymes))
    graph, enzymes_sites, core_enzymes = \
        _enzymes_names_to_distances_graph(enzymes_names)
    def one_try():
        pickable_enzymes_not_in_seq = set(enzymes_names)
        enzymes_not_in_seq = set(enzymes_names)
        enzymes_in_seq = set()
        
        def add_enzyme(e, as_core=True):
            if e not in enzymes_in_seq:
                if e in pickable_enzymes_not_in_seq:
                    pickable_enzymes_not_in_seq.remove(e)
                enzymes_not_in_seq.remove(e)
                enzymes_in_seq.add(e)
                if as_core:
                    for enz in core_enzymes.get(e, []):
                        add_enzyme(enz, as_core=False)
        
        l = list(pickable_enzymes_not_in_seq)
        path = [l[np.random.randint(len(l))]]
        seq = _enzyme_path_to_sequence(path, graph, enzymes_sites)
        add_enzyme(path[-1])
        while len(pickable_enzymes_not_in_seq):
            nodes_distances = sorted([
                (graph[(path[-1], n)]['dist'], n)
                for n in pickable_enzymes_not_in_seq
            ], key=lambda dn: (dn[0], np.random.rand()))
            choices = [n for d, n in nodes_distances]
            for new_node in choices:
                new_path = path + [new_node]
                new_seq = _enzyme_path_to_sequence(
                    new_path, graph, enzymes_sites)
                new_diff = len(new_seq) - len(seq)
                end_seq = new_seq[-(6 + new_diff):]
                matches = find_patterns_matches(end_seq, patterns)
                if unique_sites:
                    new_matches = {
                        e: [m for m in matches[e]
                            if m.end > len(end_seq) - new_diff]
                        for e in matches
                    }
                    forbidden  = list(enzymes_in_seq) + list(forbidden_enzymes)
                    if any ([len(new_matches[e]) for e in forbidden]):
                        continue
                    if any ([len(new_matches[e]) > 1
                             for e in enzymes_not_in_seq]):
                        continue
                seq = new_seq
                path = new_path
                add_enzyme(new_node)
                for e in list(enzymes_not_in_seq):
                    if len(matches[e]) > 0:
                        add_enzyme(e)
                break
            else:
                return seq, enzymes_in_seq, enzymes_not_in_seq
        return seq, enzymes_in_seq, enzymes_not_in_seq
    
    failure_score = (10000, 10000)
    def score(seq, in_seq, leftover):
        if unique_sites or len(forbidden_enzymes):
            matches = find_patterns_matches(seq, patterns)
            if unique_sites and any([len(matches[e]) != 1 for e in in_seq]):
                return failure_score
            if any([len(matches[e]) for e in forbidden_enzymes]):
                return failure_score
        if success_condition is not None:
            if not success_condition(seq, in_seq, leftover):
                return failure_score
        return (len(leftover), len(seq))
        
    current_best = seq, in_seq, leftover = one_try()
    current_score = score(seq, in_seq, leftover)
    for _ in logger.iter_bar(tries=range(tries)):
        new_result = seq, in_seq, leftover = one_try()
        new_score = score(seq, in_seq, leftover)
        if new_score < current_score:
            current_score = new_score
            current_best = new_result
    if current_score == failure_score:
        raise NoSuitableSequenceFound("Could not find a suitable sequence.")
    return current_best