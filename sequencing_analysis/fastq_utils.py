import os, sys
import math
import pickle

import numpy as np
import tqdm

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import utils

def analyze_fastq(path_forward, path_reverse, path_templates, output_prefix, max_reads = None, max_hash_defect = 5, min_cluster_abundance = 2e-4, sequencing_accuracy = 0.75, max_reconstruction_depth = 100, sample_id = None, pbar = None):
  # Load sequencing data
  reads_all = utils.import_fastq_2way(path_forward, path_reverse, max_reads = max_reads, pbar = pbar, sample_id = sample_id, max_sequence_length = None)

  # Load template data
  templates = utils.import_fasta(path_templates, pbar = pbar, sample_id = sample_id)

  # Cluster reads based on their hash sequences
  min_cluster_size = math.ceil(min_cluster_abundance * len(reads_all))
  cluster_hashes, cluster_reads = cluster_reads_by_hash(
      reads_all, 
      max_hash_defect = max_hash_defect,
      min_cluster_size = min_cluster_size,
      sample_id = sample_id, 
      pbar = pbar
  )

  # Count cluster sizes by UMI
  cluster_sizes_UMI = calc_cluster_sizes_UMI(reads_all, cluster_reads)

  # Attempt sequence reconstruction of each cluster
  reconstructions, base_identity_counts, reconstruction_depths = cluster_sequence_reconstruction(
      reads_all,
      cluster_reads,
      sequencing_accuracy = sequencing_accuracy, 
      max_reconstruction_depth = max_reconstruction_depth,
      sample_id = sample_id,
      pbar = pbar
  )

  # Attempt to associate each cluster with one of the templates
  cluster_assignments = assign_clusters_to_templates(cluster_hashes, templates, sample_id = sample_id, pbar = pbar)
  template_counts = [len(cluster_reads[cluster_assignments.index(i)]) if i in cluster_assignments else 0 for i in range(len(templates))]
 
  # Save raw analysis results
  outpath = open('{}_data.p'.format(output_prefix), 'wb')
  output = {
    'sample_id': sample_id,
    'path_forward': path_forward,
    'path_reverse': path_reverse,
    'templates': templates,
    'cluster_hashes': cluster_hashes,
    'cluster_reads': cluster_reads,
    'cluster_sizes': [len(lst) for lst in cluster_reads],
    'cluster_sizes_UMI': cluster_sizes_UMI,
    'cluster_assignments': cluster_assignments,
    'cluster_reconstructions': reconstructions,
    'cluster_reconstruction_depths': reconstruction_depths,
    'cluster_base_counts': base_identity_counts,
    'template_counts': template_counts,
    'max_reads': max_reads,
    'num_reads': len(reads_all),
    'max_hash_defect': max_hash_defect,
    'min_cluster_abundance': min_cluster_abundance,
    'sequencing_accuracy': sequencing_accuracy,
    'max_reconstruction_depth': max_reconstruction_depth,
  }
  pickle.dump(output, outpath)
  outpath.close()

  # Make plots
  if sample_id is None:
    plot_title = 'Counts from {} reads'.format(len(reads_all))
  else:
    plot_title = 'Counts from {} reads (SAMPLE {})'.format(len(reads_all), sample_id)
  plot_counts_bar(templates, template_counts, output_prefix, title=plot_title)

  return output

def cluster_reads_by_hash(reads, max_hash_defect = 5, min_cluster_size = 3, cluster_refresh_interval = 20000, sample_id = None, pbar = None):
  # Steps:
  # For each forward/reverse read, do the following:
  #  1) Identify the forward hash, which is located at forward_read[23:47]
  #  2) Identify the reverse hash, which is located at reverse_read[24:48]
  #  3) Concatenate the forward/reverse hashes to get an overall hash. Determine if this hash matches a previously observed hash
  #    a) if not, add to a list of hashes.
  #    b) if so, add the read name to a list of reads associated with that cluster
  #  4) Filter clusters that are smaller than the minimum cluster size
  # Return two lists:
  #  - a list of the hashes matched to each cluster
  #  - a list of the reads matched to each cluster

  # An issue with this approach is that some clusters may merge
  # as more sequences are processed. Thus, we periodically check
  # the cluster consensus hashes and combine clusters within
  # the allowed hash defect. We also clean up the cluster consensus
  # sequences at this time.

  # After an initial pass through the dataset to determine the clusters,
  # the dataset is processed again with fixed clusters to perform
  # a final assignment of reads to each cluster.

  # We maintain the following information throughout
  #  - a dict mapping hash subsequences to a "cluster index"
  #  - the sets of hash sequences associated with each cluster
  #  - the sets of read indices associated with each cluster
  num_clusters = 0
  hash_subseqs_to_cluster = {}
  cluster_hashes = []
  cluster_reads = []
  cluster_hash_base_counts = []
  cluster_sequences = []

  hash_subseq_length = 48//(max_hash_defect+1)
  
  for read_idx, (_, (seq_f, seq_r)) in enumerate(reads):
    # extract hash; hash_f/hash_r starts at 23/24, but UMIs add 12nts to each end
    hash_full = read_to_hash(seq_f, seq_r, hash_f_start = 35, hash_r_start = 36)
    if len(hash_full) != 48:  continue

    # check if we've observed this hash before, or one with up to 3 mutations
    potential_clusters = set(
        idx
        for i in range(max_hash_defect+1)
        for idx in hash_subseqs_to_cluster.get(hash_full[hash_subseq_length*i:hash_subseq_length*(i+1)], [])
    )
    clusters = set()
    for idx in potential_clusters:
      seq = cluster_sequences[idx]
      if utils.hamming_distance(seq, hash_full) < max_hash_defect:
        clusters.add(idx)

    # no cluster was matched, make a new (blank) cluster
    if len(clusters) == 0:
      clusters = set([num_clusters])
      cluster_hashes.append(set())
      cluster_reads.append([])
      cluster_hash_base_counts.append([[0,0,0,0] for _ in range(48)])
      cluster_sequences.append('N'*48)
      num_clusters += 1

    # associate each substring of the hash with each matching cluster
    for offset in range(0, len(hash_full) - hash_subseq_length + 1, hash_subseq_length):
      hash_subseq = hash_full[offset:offset+hash_subseq_length]
      hash_subseqs_to_cluster[hash_subseq] = hash_subseqs_to_cluster.get(hash_subseq, set()) | clusters

    # update cluster information for each matching cluster
    for idx in clusters:
      cluster_hashes[idx].add(hash_full)
      cluster_reads[idx].append(read_idx)
      for i,nt in enumerate(hash_full):
        if nt in 'ATCG':
          cluster_hash_base_counts[idx][i]['ATCG'.index(nt)] += 1
      cluster_sequences[idx] = ''.join(['ATCG'[counts.index(max(counts))] for counts in cluster_hash_base_counts[idx]])

    # If we've reached the cluster refresh interval, merge similar clusters
    if read_idx % cluster_refresh_interval == cluster_refresh_interval-1 or read_idx == len(reads)-1:
      # remap cluster indices to a single cluster index if the consensus sequences are similar enough
      new_clusters = {}
      num_new_clusters = 0
      for idx, seq in enumerate(cluster_sequences):
        matches = [i for i in range(idx) if utils.hamming_distance(cluster_sequences[i],seq) < max_hash_defect]
        if len(matches) == 0:
          new_clusters[idx] = num_new_clusters
          num_new_clusters += 1
        else:
          for i in matches + [idx]:
            new_clusters[i] = new_clusters[matches[0]]

      # remake/update cluster info
      cluster_reads_new = [[] for _ in range(num_new_clusters)]
      cluster_hash_base_counts_new = [[[0,0,0,0] for __ in range(48)] for _ in range(num_new_clusters)]
      for c_idx in new_clusters:
        c_idx_new = new_clusters[c_idx]
        cluster_reads_new[c_idx_new].extend(cluster_reads[c_idx])
        for i,(a,t,c,g) in enumerate(cluster_hash_base_counts[c_idx]):
          atot, ttot, ctot, gtot = cluster_hash_base_counts_new[c_idx_new][i]
          cluster_hash_base_counts_new[c_idx_new][i] = [a+atot, t+ttot, c+ctot, g+gtot]
      cluster_reads = [list(set(r)) for r in cluster_reads_new]
      cluster_hash_base_counts = cluster_hash_base_counts_new

      cluster_sequences = [''.join(['ATCG'[counts.index(max(counts))] for counts in base_counts]) for base_counts in cluster_hash_base_counts]

      cluster_hashes = [set([seq]) for seq in cluster_sequences]
      hash_subseqs_to_cluster = {}
      for c_idx,h in enumerate(cluster_sequences):
        for offset in range(0, len(h) - hash_subseq_length + 1, hash_subseq_length):
          h_sub = h[offset:offset+hash_subseq_length]
          hash_subseqs_to_cluster[h_sub] = hash_subseqs_to_cluster.get(h_sub, set()) | set([c_idx])

      num_clusters = num_new_clusters
          
  print("{}: Clustering reads complete.".format(sample_id))
      
  # calculate consensus hashes for each cluster
  cluster_hash_seqs = []
  for cluster_idx in range(len(cluster_reads)):
    cluster_hash_seqs.append(''.join(['ATCG'[counts.index(max(counts))] for counts in cluster_hash_base_counts[cluster_idx]]))


  # with these clusters fixed, assign reads to clusters
  cluster_reads2 = [[] for _ in range(len(cluster_hash_seqs))]
  for read_idx, (read_name, (seq_f, seq_r)) in enumerate(reads):
    # extract hash; hash_f/hash_r starts at 23/24, but UMIs add 12nts to each end
    hash_full = read_to_hash(seq_f, seq_r, hash_f_start = 35, hash_r_start = 36)
    if len(hash_full) != 48:  continue

    dists = [utils.hamming_distance(hash_full, h) for h in cluster_hash_seqs]
    min_dist = min(dists)
    if min_dist <= max_hash_defect:
      cluster_reads2[dists.index(min_dist)].append(read_idx)

  # filter clusters that are too small
  cluster_hash_seqs_filt, cluster_reads_filt = [],[]
  for hash_seq, read_idxs in zip(cluster_hash_seqs, cluster_reads2):
    if len(read_idxs) >= min_cluster_size:
      cluster_hash_seqs_filt.append(hash_seq)
      cluster_reads_filt.append(read_idxs)
  
  print("{}: Assignment of reads to clusters complete.".format(sample_id))
 
  if len(cluster_hash_seqs_filt) > 0:
    return zip(*sorted(zip(cluster_hash_seqs_filt, cluster_reads_filt), key=lambda x: len(x[1]), reverse=True))
  else:
    return [],[]

def calc_cluster_sizes_UMI(reads, cluster_read_indices):
  # conflates reads that have UMIs that are different by up to 1 mutation

  cluster_sizes_UMI = []

  for read_idxs in cluster_read_indices:
    if len(read_idxs) == 0:
      cluster_sizes_UMI.append(0)
      continue

    UMI_counts = {}
    for read_idx in read_idxs:
      _, (seq_f, seq_r) = reads[read_idx]
      UMI = read_to_UMI(seq_f, seq_r)
      UMI_counts[UMI] = UMI_counts.get(UMI, 0) + 1

    UMIs_sorted,_ = zip(*sorted(UMI_counts.items(), key = lambda x: x[1], reverse=True))

    UMIs_filt = set()
    UMIs_used = set()
    for UMI in UMIs_sorted:
      if UMI not in UMIs_used:
        UMIs_filt.add(UMI)
        UMIs_used.add(UMI)
        for i in range(len(UMI)):
          for nt in 'ATCG':
            UMIs_used.add(UMI[:i] + nt + UMI[i+1:])

    cluster_sizes_UMI.append(len(UMIs_filt))

  return cluster_sizes_UMI


def cluster_sequence_reconstruction(reads, cluster_read_indices, sequencing_accuracy = 0.75, max_reconstruction_depth = 100, sample_id = None, pbar = None):

  # Steps:
  # For each cluster:
  #  1) For each read in this cluster
  #    a) Pull out the corresponding forward and reverse sequences
  #    b) Attempt to align the forward and reverse sequences
  #    c) If sequence alignment is successful, update:
  #      - the number of aligned sequences used in the reconstruction
  #      - the base identity counts at each nt position
  #    d) Stop processing reads early if more than <max_reconstruction_depth> sequences have aligned successfully
  #  2) Perform a majority-vote consensus base at each nucleotide location
  #    - the votes between ATCG and no base are compared; if the "no base" wins out, the sequence
  #      is assumed to end earlier than the longest read used in the alignment

  reconstruction_depths = []
  base_identity_counts_all = []
  consensus_sequences = []
  for read_idxs in cluster_read_indices:

    reconstruction_depth = 0
    base_identity_counts = []
    for read_idx in read_idxs:
      _, (seq_f, seq_r) = reads[read_idx]
     
      # attempt to align the forward and reverse reads
      # this is done by taking the last 20nts of the reverse read and aligning it with the forward read
      seq_r_sub = utils.sequence_complement(seq_r[-20:])
      align_pos, align_score = utils.align_sequences(seq_r_sub, seq_f)
      if align_score >= len(seq_r_sub)*sequencing_accuracy:
        # alignment successful
        reconstruction_depth += 1

        seq_full = seq_f + utils.sequence_complement(seq_r)[len(seq_f)-align_pos:]
        if len(base_identity_counts) < len(seq_full):
          base_identity_counts += [(0,0,0,0)]*(len(seq_full) - len(base_identity_counts))
        for i,nt in enumerate(seq_full):
          a,t,c,g = base_identity_counts[i]
          if   nt=='A':  a+=1
          elif nt == 'T':  t+=1
          elif nt == 'C':  c+=1
          elif nt == 'G':  g+=1
          base_identity_counts[i] = (a,t,c,g)

      if reconstruction_depth >= max_reconstruction_depth:
        break

    consensus_sequence = ''
    for counts in base_identity_counts:
      no_base = reconstruction_depth - sum(counts)
      max_nt = max(counts)

      if no_base > max_nt:  break
      
      if counts.count(max_nt) > 1:
        consensus_sequence += 'N'
      else:
        consensus_sequence += 'ATCG'[counts.index(max_nt)]

    reconstruction_depths.append(reconstruction_depth)
    base_identity_counts_all.append(base_identity_counts)
    consensus_sequences.append(consensus_sequence)

  print("{}: Sequence reconstruction complete.".format(sample_id))

  return consensus_sequences, base_identity_counts_all, reconstruction_depths

def assign_clusters_to_templates(cluster_hashes, templates, max_hash_defect = 3, sample_id = None, pbar = None):
  # Extract hashes for each template:
  template_hashes = []
  for _, template_seq in templates:
    template_hashes.append(sequence_to_hash(template_seq, hash_f_start = 23, hash_r_start = 24))

  # For each cluster, see if one of the template hashes is included in the cluster's list of hashes
  cluster_assignments = []
  for h in cluster_hashes:
    dists = [utils.hamming_distance(t_hash, h) for t_hash in template_hashes]
    min_dist = min(dists)
    if min_dist <= max_hash_defect:
      cluster_assignments.append(dists.index(min_dist))
    else:
      cluster_assignments.append(None)

  return cluster_assignments

def plot_counts_bar(templates, counts, output_prefix, title = None):

  bar_x = np.arange(len(templates))
  bar_heights = counts
  bar_width = .9
  bar_tick_labels = [name[:name.index('_')] for name,_ in templates]
  ylims = (0, max(counts)*1.02+.01)

  plt.figure()
  plt.bar(x = bar_x, height = bar_heights, width = bar_width, data = counts)
  plt.ylabel('Counts')
  plt.xticks(bar_x, bar_tick_labels, rotation='vertical')

  for x, c in zip(bar_x, counts):
    plt.text(x, ylims[1]/50., str(c), horizontalalignment = 'center', verticalalignment = 'bottom', rotation='vertical')

  plt.ylim(ylims)

  if title is not None:
    plt.title(title)

  plt.savefig('{}_counts.pdf'.format(output_prefix), bbox_inches='tight')


def read_to_hash(seq_f, seq_r, hash_f_start = 23, hash_f_length = 24, hash_r_start = 24, hash_r_length = 24):
  # Forward hash is given by:
  #  seq_f[hash_f_start:hash_f_start+hash_f_length] 
  # Reverse hash is given by:
  #  seq_r[hash_r_start:hash_r_start+hash_r_length]
  # Full hash is:
  #  forward_hash + complement(reverse_hash)
  hash_f = seq_f[hash_f_start:hash_f_start+hash_f_length]
  hash_r = seq_r[hash_r_start:hash_r_start+hash_r_length]
  return hash_f + utils.sequence_complement(hash_r)

def sequence_to_hash(sequence, hash_f_start = 23, hash_f_length = 24, hash_r_start = 24, hash_r_length = 24):
  # Parameters are the same as for read_to_hash
  seq_r = utils.sequence_complement(sequence)
  return read_to_hash(sequence, seq_r, hash_f_start, hash_f_length, hash_r_start, hash_r_length)

def read_to_UMI(seq_f, seq_r, UMI_f_start = 0, UMI_f_length = 12, UMI_r_start = 0, UMI_r_length = 12):
  # Parameters are interpreted as in read_to_hash()
  # but default UMI lengths are 12 for both forward and reverse
  UMI_f = seq_f[UMI_f_start:UMI_f_start+UMI_f_length]
  UMI_r = seq_r[UMI_r_start:UMI_r_start+UMI_r_length]
  return UMI_f + utils.sequence_complement(UMI_r)


