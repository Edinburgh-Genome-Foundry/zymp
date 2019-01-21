from copy import deepcopy
from Bio import Restriction, SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dnachisel import DnaNotationPattern
import proglog

def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses BioPython for speed.
    """
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.

    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]

def annotate_record(seqrecord, location="full", feature_type="misc_feature",
                    margin=0, **qualifiers):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord)-margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type
        )
    )

def sequence_to_biopython_record(sequence, id='<unknown id>',
                                 name='<unknown name>', features=()):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(Seq(sequence, alphabet=DNAAlphabet()),
                     id=id, name=name, features=list(features))

def enzymes_to_dna_pattern(enzymes_names):
    return {
        e: DnaNotationPattern(Restriction.__dict__[e].site,
                              in_both_strands='auto')
        for e in enzymes_names
    }

def get_enzymes_ATGC_sequences(enzymes_names):
    patterns = enzymes_to_dna_pattern(enzymes_names)
    return {e: p.all_variants()[0] for e, p in patterns.items()}

def find_patterns_matches(seq, patterns=None, enzymes_names=None):
    if enzymes_names is not None:
        patterns = enzymes_to_dna_pattern(enzymes_names)
    return {
        e: pattern.find_matches(seq)
        for e, pattern in patterns.items()
    }

def annotate_enzymes_sites(sequence, enzymes_names, forbidden_enzymes=(),
                           unique_sites=True, valid_color='#ccccff',
                           invalid_color='#ffcccc'):
    record = sequence_to_biopython_record(sequence)
    matches = find_patterns_matches(
        sequence.upper(), enzymes_names=enzymes_names)
    for enz_name, site_matches in matches.items():
        if unique_sites and len(site_matches) > 1:
            color = invalid_color
        else:
            color = valid_color
        for m in site_matches:
            annotate_record(record, m.to_tuple(), label=enz_name, color=color)
    if len(forbidden_enzymes):
        matches = find_patterns_matches(sequence,
                                        enzymes_names=forbidden_enzymes)
        for enz_name, site_matches in matches.items():
            for m in site_matches:
                annotate_record(record, m.to_tuple(), label=enz_name,
                                color='red')
    return record

def write_record(record, target, fmt='genbank'):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes"""
    record = deepcopy(record)
    record.name = record.name[:20]
    if str(record.seq.alphabet.__class__.__name__) != 'DNAAlphabet':
        record.seq.alphabet = DNAAlphabet()
    if hasattr(target, 'open'):
        target = target.open('w')
    SeqIO.write(record, target, fmt)
