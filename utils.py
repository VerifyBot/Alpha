import json
import os
from itertools import chain
import Bio.Data.CodonTable as ct
from scipy.stats import gmean
from collections import Counter

# get rid of Biopython warning
import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)


def count_nucs_sliding(rand_nucs: str, seq: str = None, nucs_count: int = 3, frame: int = 0, on='variant') -> dict:
  """
  Count how many `nucs_count` length nucs and which appear
  in the `seq` sequence in a sliding mode
  """
  counter = {on: seq or rand_nucs}

  rand_nucs = rand_nucs[frame:]

  for i in range(0, len(rand_nucs) - nucs_count + 1):
    nucs = rand_nucs[i:i + nucs_count]
    cn = nucs + '_count' + str(frame if frame != 0 else '')

    if cn in counter:
      counter[cn] += 1
    else:
      counter[cn] = 1

  return counter


def count_codons(sequence, is_sliding=True, codon_length=3):
  """
  Count how many of each codon appear in the sequence either
  in a sliding mode or in a slicing mode (is_sliding=False)
  """
  assert len(sequence) > 3

  d = {}
  i = 0
  while i <= len(sequence) - codon_length:
    codon = sequence[i:i + codon_length]
    if codon in d:
      d[codon] += 1
    else:
      d[codon] = 1
    if is_sliding:
      i += 1
    else:
      i += 3
  return d


def apply_counts(row):
  """
  Function used to apply func::count_codons on the variants dataframe
  for each of its parts (variant, prefix, rand_nucs, suffix) in several frames
  and also in a sliding mode
  """

  for cn in ['variant', 'prefix', 'rand_nucs', 'suffix']:
    for method in ['sliding', 'slicing_frame0', 'slicing_frame1', 'slicing_frame2']:
      is_slicing = method.startswith('slicing')
      string_to_work_on = row[cn]
      if is_slicing:
        string_to_work_on = string_to_work_on[int(method[-1]):]
      frequency_counts_dictionary = count_codons(string_to_work_on, is_slicing)
      for k, v in frequency_counts_dictionary.items():
        row[f"{cn}_{method}_{k}_count"] = v
  return row


BASE_FOLDER = 'Alpha'
fpath = '../' if not os.getcwd().endswith(BASE_FOLDER) else ''
with open(fpath + 'data/synonymous_codons.json') as f:
  _synonymous_codons = json.load(f)
_non_synonymous_codons = {'ATG', 'TGG'}


def RSCU(sequences, genetic_code=11):
  r"""Calculates the relative synonymous codon usage (RSCU) for a set of sequences.
  RSCU is 'the observed frequency of [a] codon divided by the frequency
  expected under the assumption of equal usage of the synonymous codons for an
  amino acid' (page 1283).
  In math terms, it is
  .. math::
      \frac{X_{ij}}{\frac{1}{n_i}\sum_{j=1}^{n_i}x_{ij}}
  "where :math:`X` is the number of occurrences of the :math:`j` th codon for
  the :math:`i` th amino acid, and :math:`n` is the number (from one to six)
  of alternative codons for the :math:`i` th amino acid" (page 1283).
  Args:
      sequences (list): The reference set of sequences.
      genetic_code (int, optional): The translation table to use. Defaults to 11, the standard genetic code.
  Returns:
      dict: The relative synonymous codon usage.
  Raises:
      ValueError: When an invalid sequence is provided or a list is not provided.
  """

  if not isinstance(sequences, (list, tuple)):
    raise ValueError(
      "Be sure to pass a list of sequences, not a single sequence. "
      "To find the RSCU of a single sequence, pass it as a one element list."
    )

  # ensure all input sequences are divisible by three
  for sequence in sequences:
    if len(sequence) % 3 != 0:
      raise ValueError("Input sequence not divisible by three")
    if not sequence:
      raise ValueError("Input sequence cannot be empty")

  # count the number of each codon in the sequences
  sequences = (
    (sequence[i: i + 3].upper() for i in range(0, len(sequence), 3))
    for sequence in sequences
  )
  codons = chain.from_iterable(
    sequences
  )  # flat list of all codons (to be used for counting)
  counts = Counter(codons)

  # "if a certain codon is never used in the reference set... assign [its
  # count] a value of 0.5" (page 1285)
  for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
    if counts[codon] == 0:
      counts[codon] = 0.5

  # determine the synonymous codons for the genetic code
  synonymous_codons = _synonymous_codons

  # hold the result as it is being calulated
  result = {}

  # calculate RSCU values
  for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
    result[codon] = counts[codon] / (
        (len(synonymous_codons[codon]) ** -1)
        * (sum((counts[_codon] for _codon in synonymous_codons[codon])))
    )

  return result


def relative_adaptiveness(sequences=None, RSCUs=None, genetic_code=11):
  r"""Calculates the relative adaptiveness/weight of codons.
  The relative adaptiveness is "the frequency of use of that codon compared to
  the frequency of the optimal codon for that amino acid" (page 1283).
  In math terms, :math:`w_{ij}`, the weight for the :math:`j` th codon for
  the :math:`i` th amino acid is
  .. math::
      w_{ij} = \frac{\text{RSCU}_{ij}}{\text{RSCU}_{imax}}
  where ":math:`\text{RSCU}_{imax}` [is] the RSCU... for the frequently used
  codon for the :math:`i` th amino acid" (page 1283).
  Args:
      sequences (list, optional): The reference set of sequences.
      RSCUs (dict, optional): The RSCU of the reference set.
      genentic_code (int, optional): The translation table to use. Defaults to 11, the standard genetic code.
  Note:
      Either ``sequences`` or ``RSCUs`` is required.
  Returns:
      dict: A mapping between each codon and its weight/relative adaptiveness.
  Raises:
      ValueError: When neither ``sequences`` nor ``RSCUs`` is provided.
      ValueError: See :func:`RSCU` for details.
  """

  # ensure user gave only and only one input
  if sum([bool(sequences), bool(RSCUs)]) != 1:
    raise TypeError("Must provide either reference sequences or RSCU dictionary")

  # calculate the RSCUs if only given sequences
  if sequences:
    RSCUs = RSCU(sequences, genetic_code=genetic_code)

  # determine the synonymous codons for the genetic code
  synonymous_codons = _synonymous_codons

  # calculate the weights
  weights = {}
  for codon in RSCUs:
    weights[codon] = RSCUs[codon] / max(
      (RSCUs[_codon] for _codon in synonymous_codons[codon])
    )

  return weights


def CAI(sequence, weights):
  """CAI Index calculation"""
  sequence = sequence.upper()
  sequence = [sequence[i: i + 3] for i in range(0, len(sequence), 3)]

  sequence_weights = []
  for codon in sequence:
    if codon not in _non_synonymous_codons:
      try:
        sequence_weights.append(weights[codon])
      except KeyError:
        # ignore stop codons
        if codon in ct.unambiguous_dna_by_id[11].stop_codons:
          pass

  # return the geometric mean of the weights raised to one over the length of the sequence
  return float(gmean(sequence_weights))


def TAI(sequence, weights):
  """TAI Index is the same as CAI calculation except for end codon being irrelevant"""
  return CAI(sequence[:-3], weights)


codons_amino_acid_table = {
  'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
  'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
  'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

# Amino acid degeneracy
degeneracy = {2: ['M'], 9: ['N', 'D', 'C', 'Q', 'E', 'H', 'K', 'F', 'Y'], 1: ['I'], 5: ['A', 'G', 'P', 'T', 'V'],
              3: ['S', 'L', 'R']}


def get_deg(aa):
  for k, v in degeneracy.items():
    if aa in v:
      return k


def ENC(sequence: str) -> int:
  """ENC Index -- /shrug"""
  amino_acids = {}

  for i in range(0, len(sequence), 3):
    codon = sequence[i:i + 3]
    amino_acid = codons_amino_acid_table[codon]

    # appearTimes += 1
    if amino_acid in amino_acids:
      amino_acids[amino_acid]['appearTimes'] += 1
    else:
      amino_acids[amino_acid] = {'codons': {}, 'appearTimes': 1}

    # codon in aa += 1
    if codon in amino_acids[amino_acid]['codons']:
      amino_acids[amino_acid]['codons'][codon] += 1
    else:
      amino_acids[amino_acid]['codons'][codon] = 1

  # calculate Pi
  for aa, amino_data in amino_acids.items():
    for codon in amino_data['codons'].keys():
      amino_acids[aa]['codons'][codon] /= amino_acids[aa]['appearTimes']

  # calculate F
  for aa, amino_data in amino_acids.items():
    amino_acids[aa]['F'] = sum(v ** 2 for v in amino_acids[aa]['codons'].values())

  # calculate ENC
  return sum(get_deg(k) / (1 / v['F']) for k, v in amino_acids.items() if k not in ['W', '_'])
