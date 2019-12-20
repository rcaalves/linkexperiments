import re, sys
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.stats


def parsefile(filename):
  f = file(filename)

  # Line example: 1576792544 E5.3B>C7.49| 304L106R73
  expression = re.compile(r"""^
    (?P<TIME>\d+)\s
    (?P<SRC>[0-9a-fA-F]{2}\.[0-9a-fA-F]{2})>
    (?P<DST>[0-9a-fA-F]{2}\.[0-9a-fA-F]{2})\|
    \s*(?P<SEQNO>\d+)
    L(?P<LQI>\d+)
    R(?P<RSSI>\d+)$
    """, re.X) #re.X: verbose so we can comment along

  ts_sequence = []
  rssi_sequence = []
  lqi_sequence = []
  seqno_sequence = []

  previous_src = previous_dst = ""

  button_press = False

  for line in f:
    s = expression.search(line)
    timestamp, src, dst, seqno, lqi, rssi = (None,) * 6

    if line.strip() == "BTN":
      button_press = True

    if s != None:
      # print line.strip()
      timestamp, src, dst, seqno, lqi, rssi = s.groups()
      seqno, lqi, rssi = map(int, (seqno, lqi, rssi))

    # print timestamp, src, dst, seqno, lqi, rssi
    # print

      ts_sequence += [timestamp]
      rssi_sequence += [rssi]
      lqi_sequence += [lqi]
      seqno_sequence += [seqno]

      if previous_src != src and previous_src not in (None, "") and src not in (None, ""):
        print >> sys.stderr, "Unexpected change of source address", previous_src, ">", src

      if previous_dst != dst and previous_dst not in (None, "") and dst not in (None, ""):
        print >> sys.stderr, "Unexpected change of destination address", previous_dst, ">", dst

      previous_src = src
      previous_dst = dst
  f.close()

  return seqno_sequence, lqi_sequence, rssi_sequence, button_press

def check():
  seqno_sequence = filter(None, seqno_sequence)
  print seqno_sequence
  i = 1
  while i < len(seqno_sequence):
    if seqno_sequence[i] > seqno_sequence[i - 1] + 1:
      missed_packets = seqno_sequence[i] - seqno_sequence[i - 1] - 1
      print "miss", seqno_sequence[i-1], missed_packets

    i+= 1

def fill_with_none_and_cap(seq):
  # i = 1
  # seq_temp = [seq[0]]
  # while i < len(seq):
  #   if seq[i] < seq[i-1]:
  #     break
  #
  #   seq_temp += [None] * (seq[i] - seq[i-1] - 1)
  #   seq_temp += [seq[i]]
  #
  #   i+= 1

  prev = seq[0]
  seq_temp = [prev]
  for i in seq[1:]:
    if i < prev:
      break

    seq_temp += [None] * (i - prev - 1)
    seq_temp += [i]

    prev = i
  return seq_temp

def insert_nones(addy, model):
  i = 0
  while i < len(model):
    if model[i] == None and addy[i] != None:
      addy.insert(i, None)
    i+=1
  while len(addy) > len(model):
    addy.pop()

#
def pair_up(seqs0, seqs1):
  if seqs1[-1] == True:
    seqs0, seqs1 = seqs1, seqs0

  s0, l0, r0 = seqs0[:3]
  s1, l1, r1 = seqs1[:3]

  s0 = fill_with_none_and_cap(s0)
  s1 = fill_with_none_and_cap(s1)
  s0 = [None] * (s0[0] - 1) + s0

  insert_nones(l0, s0)
  insert_nones(l1, s1)

  insert_nones(r0, s0)
  insert_nones(r1, s1)

  # print zip (s0, s1)[-10:], len(s0), len(s1)
  # print l0[-10:], len(l0)
  # print l1[-10:], len(l1)

  return zip(s0, s1), zip(l0, l1), zip(r0, r1)

n0 = parsefile("node0.log")
n1 = parsefile("node1.log")

seq_no_pairs, lqi_pairs, rssi_pairs = pair_up(n0, n1)
plt.scatter(*zip(*lqi_pairs))
# plt.show()
lqi_pairs = filter(lambda x: None not in x, lqi_pairs)
lqi_coef = np.corrcoef([x[0] for x in lqi_pairs], [x[1] for x in lqi_pairs])[0][1]
print "pearson coef LQI", lqi_coef

plt.scatter(*zip(*rssi_pairs))
# plt.show()
rssi_pairs = filter(lambda x: None not in x, rssi_pairs)
rssi_coef = np.corrcoef([x[0] for x in rssi_pairs], [x[1] for x in rssi_pairs])[0][1]
print "pearson coef RSSI", rssi_coef
