# compatilibity modules
from __future__ import print_function
from future.utils import lfilter, lzip, lrange, lmap
from io import open

import re, sys
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy.stats
import random

def parsefile(filename):
  f = open(filename, 'r', encoding="latin1")

  # Line example: 1576792544 E5.3B>C7.49| 304L106R73
  expression = re.compile(r"""^
    (?P<TIME>\d+)\s
    (?P<SRC>[0-9a-fA-F]{1,2}\.[0-9a-fA-F]{1,2})>
    (?P<DST>[0-9a-fA-F]{1,2}\.[0-9a-fA-F]{1,2})\|
    \s*(?P<SEQNO>\d+)
    L(?P<LQI>\d+)
    R(?P<RSSI>\d+)$
    """, re.X) #re.X: verbose, so we can comment along

  ts_sequence = []
  rssi_sequence = []
  lqi_sequence = []
  seqno_sequence = []

  previous_src = previous_dst = ""

  extra_info = []

  previous_seqno = None

  for line in f:
    s = expression.search(line)
    timestamp, src, dst, seqno, lqi, rssi = (None,) * 6
    line_no_ts = line.strip()[line.strip().find(" "):]

    if "BTN" in line_no_ts:
      if "BTN" in extra_info:
        break
      extra_info += ["BTN"]

    if "INIT" in line_no_ts:
      init_number = int(line_no_ts[line_no_ts.find("INIT ")+len("INIT "):])+1
      extra_info += [init_number, previous_seqno if previous_seqno != None else 0]

    if s != None:
      # print (line.strip())
      timestamp, src, dst, seqno, lqi, rssi = s.groups()
      seqno, lqi, rssi = map(int, (seqno, lqi, rssi))
      previous_seqno = seqno

    # print (timestamp, src, dst, seqno, lqi, rssi)
    # print ()

      ts_sequence += [timestamp]
      rssi_sequence += [rssi]
      lqi_sequence += [lqi]
      seqno_sequence += [seqno]

      if previous_src != src and previous_src not in (None, "") and src not in (None, ""):
        print("Unexpected change of source address", previous_src, ">", src, file=sys.stderr)

      if previous_dst != dst and previous_dst not in (None, "") and dst not in (None, ""):
        print("Unexpected change of destination address", previous_dst, ">", dst, file=sys.stderr)

      previous_src = src
      previous_dst = dst
  f.close()

  return seqno_sequence, lqi_sequence, rssi_sequence, extra_info

def average_loss_streak(s):
  loss_streaks = []
  streak = 0
  for i in s:
    if i == None:
      streak += 1
    elif streak != 0:
      loss_streaks += [streak]
      streak = 0
  if streak != 0:
    loss_streaks += [streak]
  if loss_streaks == []:
    loss_streaks = [0]
  # print (np.mean(loss_streaks), np.std(loss_streaks), loss_streaks)
  return np.mean(loss_streaks), np.std(loss_streaks)

def average_delivery_streak(s):
  delivery_streaks = []
  streak = 0
  for i in s:
    if i != None:
      streak += 1
    elif streak != 0:
      delivery_streaks += [streak]
      streak = 0
  if streak != 0:
    delivery_streaks += [streak]
  if delivery_streaks == []:
    delivery_streaks = [0]
  # print (np.mean(delivery_streaks), np.std(delivery_streaks), delivery_streaks)
  return np.mean(delivery_streaks), np.std(delivery_streaks)

kpe_dict_values = {}
MAX_KPE_DICT_MEM = 8*1024*1024*1024

def kpe_dict2(t, k_plus_vector, k_minus_vector):
  k_plus_t = 0
  k_minus_t = 0

  if (sys.getsizeof(kpe_dict_values) > MAX_KPE_DICT_MEM):
    kpe_dict_values.clear()

  for (m, s) in k_plus_vector:
    m_shifted, t_shifted = 0, abs(m - t)
    if (m_shifted,s,t_shifted) in kpe_dict_values:
      # print "kaching 2"
      v = kpe_dict_values[(m_shifted, s, t_shifted)]
    else:
      v = scipy.stats.norm(m_shifted, s).pdf(t_shifted)
      kpe_dict_values[(m_shifted, s, t_shifted)] = v
    k_plus_t += v

  for (m, s) in k_minus_vector:
    m_shifted, t_shifted = 0, abs(m - t)
    if (m_shifted,s,t_shifted) in kpe_dict_values:
      # print "kaching 3"
      v = kpe_dict_values[(m_shifted, s, t_shifted)]
    else:
      v = scipy.stats.norm(m_shifted, s).pdf(t_shifted)
      kpe_dict_values[(m_shifted, s, t_shifted)] = v
    k_minus_t += v

  kpe_t = k_plus_t/(k_plus_t + k_minus_t)
  return kpe_t

def kpe(t, k_plus_vector, k_minus_vector):
  return kpe_dict2(t, k_plus_vector, k_minus_vector)

double_ues = {}
def asimmetry_metric_kpe(s0, s1):
  should_print = True
  max_i = min(len(s0),len(s1))
  s0 = s0[:max_i]
  s1 = s1[:max_i]

  if max_i > 500:
    should_print = False

  s0_k_plus_gaussians = []
  s0_k_minus_gaussians = []
  s1_k_plus_gaussians = []
  s1_k_minus_gaussians = []

  if should_print:
    limcol = 20
    for j in range(0, len(s0), limcol):
      print ("   ", end="")
      for i in s0[j:min(j+limcol, len(s0))]:
        out = "X" if i == None else repr(i)
        print ("%5s" % (out,), end="")
      print ()
      print (" "*6, end="")
      for i in s1[j:min(j+limcol, len(s1))]:
        out = "X" if i == None else repr(i)
        print ("%5s" % (out,), end="")
      print ("\n")

  def cost_function(w, N):
    if w <= 0:
      return np.inf
    return N/w + (N*(N-1))/w*( np.exp(-2/(4*w*w)) - 2*math.sqrt(2)*np.exp(-2/(2*w*w)) )

  # print (scipy.optimize.minimize(cost_function, x0=1, args=(max_i,))#, method="Newton-CG"))
  # print ()
  # print (scipy.optimize.minimize_scalar(cost_function, args=(max_i,), method="Golden"))

  if max_i in double_ues:
    w = double_ues[max_i]
  else:
    min_w_result = scipy.optimize.minimize_scalar(cost_function, args=(max_i,), method="Golden")
    w = min_w_result.x
    double_ues[max_i] = w

  t = 0
  for i in s0:
    if i == None:
      s0_k_minus_gaussians += [ (t, w) ]
    else:
      s0_k_plus_gaussians += [ (t, w) ]
    t += 2

  t = 1
  for i in s1:
    if i == None:
      s1_k_minus_gaussians += [ (t, w) ]
    else:
      s1_k_plus_gaussians += [ (t, w) ]
    t += 2

  def integrand(t):
    return abs(kpe(t, s0_k_plus_gaussians, s0_k_minus_gaussians) - kpe(t, s1_k_plus_gaussians, s1_k_minus_gaussians))


  x_values = np.arange(0, max_i*2-1, 0.125)
  y_values = lmap(integrand, x_values)
  # plt.plot(x_values, y_values, marker='o', color = "black")
  # x_values2 = [x_values[i] for i in lrange(0, len(x_values), 2)]
  # y_values2 = [y_values[i] for i in lrange(0, len(y_values), 2)]
  # plt.plot(x_values2 , y_values2, marker='o', color = "green")
  # x_values2 = [x_values[i] for i in lrange(0, len(x_values), 4)]
  # y_values2 = [y_values[i] for i in lrange(0, len(y_values), 4)]
  # plt.plot(x_values2, y_values2, marker='o', color = "red")
  # x_values2 = [x_values[i] for i in lrange(0, len(x_values), 8)]
  # y_values2 = [y_values[i] for i in lrange(0, len(y_values), 8)]
  # plt.plot(x_values2, y_values2, marker='o', color = "blue")
  # plt.show()

  # max_err = 1e-2
  # m_old = scipy.integrate.quad(integrand, 0, max_i*2-1,  epsabs=max_err, epsrel=max_err)[0]/(max_i*2-1)
  m = scipy.integrate.simps(y_values, x_values)/(max_i*2-1)
  # , limit=100
  if should_print:
    print ("m kpe_dict2", m)

  if should_print:
    print ()
  print (sys.getsizeof(kpe_dict_values))
  return m

def segment_diff_integral(p1, p2, p3, p4):
  # line equation mx+c
  # m = delta(y)/delta(x)
  # c = y-mx, for a give (x,y)
  # print()
  # print(p1, p2, p3, p4)
  if (p1[0] != p3[0] or p2[0] != p4[0]):
    raise Exception("Different intervals")

  if (p2[0] <= p1[0]):
    return 0.0

  m0 = (p2[1] - p1[1])/(p2[0] - p1[0])
  c0 = p1[1]-m0*p1[0]

  m1 = (p4[1] - p3[1])/(p4[0] - p3[0])
  c1 = p3[1]-m1*p3[0]

  # print(m0,m1)
  # print(c0,c1)

  if (p1[1] > p3[1] and p2[1] < p4[1]) or (p1[1] < p3[1] and p2[1] > p4[1]):
    # print("lines are crossing")
    p_cross = ((c1-c0)/(m0-m1), (m0*c1-m1*c0)/(m0-m1))
    # print(p_cross)
    return ( segment_diff_integral(p1,      p_cross,  p3,       p_cross) +
             segment_diff_integral(p_cross, p2,       p_cross,  p4) )
  else:
    # print("lines are not crossing")
    i = p2[0]
    j = p1[0]
    value = abs( (m0-m1)/2*(i*i)+(c0-c1)*(i) -
                ((m0-m1)/2*(j*j)+(c0-c1)*(j)) )
    # print(value)
    return value


def integreate_sequence_of_segments(v0, v1):
  # v0: list of y values, corresponding to x 0, 2, 4, ...
  # v1: list of y values, corresponding to x 1, 3, 5, ...
  # returns the integral of the difference between line segments
  v0 = map(lambda x: x*1.0, v0)
  v1 = map(lambda x: x*1.0, v1)
  max_i = min(len(v0),len(v1))

  # print (v0, v1)

  acc = 0.0
  for i in lrange(1,max_i*2):

    i0 = i//2
    i1 = i//2 - 1
    # print (i, i0, i1)

    if i % 2 == 0:
      p1 = (i-1, (v0[i0] + v0[i0-1])/2)
      p2 = (i, v0[i0])
      p3 = (i-1, v1[i1])
      p4 = (i, (v1[i1] + v1[i1+1])/2)
    else:
      p1 = (i-1, v0[i0])
      p2 = (i  , (v0[i0] + v0[i0+1])/2) if i0+1 < max_i else (i, v0[i0])
      p3 = (i-1, (v1[i1] + v1[i1+1])/2) if i1 != -1 else (i-1, v1[i1+1])
      p4 = (i, v1[i1+1])

    # print (p1, p2, p3, p4)
    acc += segment_diff_integral(p1,p2,p3,p4)

    # print(acc)
    # print()

  return acc

def asimmetry_metric_simplified(s0, s1):
  max_i = min(len(s0),len(s1))
  s0 = s0[:max_i]
  s1 = s1[:max_i]

  half_window = 3
  min_weight = 0.1
  weights = []
  for j in lrange(-half_window, half_window+1):
    dist2 = j*j
    weight = (min_weight-1)/(half_window * half_window) * dist2 + 1
    # print (j, weight)
    weights += [weight]
  # div = sum(weights)
  # print (weights)

  dp0 = []
  dp1 = []
  for i in lrange(0, max_i):
    dp_temp0 = 0.0
    dp_temp1 = 0.0
    div = 0.0

    for j in lrange(max(0, i-half_window), min(max_i-1, i+half_window)+1):
      # print(j,s0[j], j-i+half_window, weights[j-i+half_window])
      div += weights[j-i+half_window]
      if s0[j] != None:
        dp_temp0 += weights[j-i+half_window]
      if s1[j] != None:
        dp_temp1 += weights[j-i+half_window]

    dp0 += [dp_temp0/div]
    dp1 += [dp_temp1/div]

  # print (map(lambda x: "%.2f" % (x,), dp0))
  # print (map(lambda x: "%.2f" % (x,), dp1))

  return integreate_sequence_of_segments(dp0, dp1)/(max_i*2)


def interpolate(d, interval=1, window=2, avg_method="Even"):
  if len(d) == 0:
    return {}, []
  max_i = max(d)

  out_dict = {}
  out_list = []

  pre_window = []
  pre_i = []

  after_i = [1]
  after_window = []
  while len(after_window) < window and after_i[-1] <= max_i:
    if after_i[-1] in d:
      after_window += [d[after_i[-1]]]
      after_i += [after_i[-1]+1]
    else:
      after_i = after_i[:-1] + [after_i[-1]+1]
  after_i = after_i[:-1]

  for i in lrange(0, (interval+1)//2 + max_i, interval):
    if i in d:
      out_dict[i] = 1.0*d[i]
      out_list += [1.0*d[i]]
    else:
      # out_dict[i] = 1.0*(sum(pre_window) + sum(after_window)) / (len(pre_window) + len(after_window))
      if avg_method == "Even":
        avg_f = lambda maxdist, dist: 1.0
      elif avg_method == "Linear":
        avg_f = lambda maxdist, dist: (1 + max_distance-abs(w-i))
      elif avg_method == "Exp":
        avg_f = lambda maxdist, dist: 2**(max_distance-abs(w-i))
      else:
        raise Exception("Unknown avg_method. Should be Even, Linear or Exp")

      out_dict[i] = 0
      div = 0
      max_distance = 1.0*max(lmap(lambda x: abs(x-i), pre_i + after_i))

      for v, w in lzip(pre_window, pre_i) + lzip(after_window, after_i):
        # w = (max_distance-abs(w-i))
        w = avg_f(max_distance, abs(w-i))
        out_dict[i] += v*w
        div += w
      out_dict[i] /= div
      out_list += [out_dict[i]]
      # calculate weighted average?

    # print (i, "=", out_dict[i], pre_i, pre_window, after_i, after_window)

    if i in d:
      if len(pre_i) == window: pre_i = pre_i[1:]
      pre_i = pre_i + [i]
      if len(pre_window) == window: pre_window = pre_window[1:]
      pre_window += [ d[pre_i[-1]] ]

    if len(after_i):
      if after_i[0] == i+1:
        after_i = after_i[1:] + [after_i[-1]+1]
        while after_i[-1] <= max_i and after_i[-1] not in d:
          after_i[-1] += 1

      if after_i[-1] not in d:
        if len(after_window):
          after_window = after_window[1:]
          after_i = after_i[:-1]
      else:
        after_window = after_window[1:] + [d[after_i[-1]]]

  return out_dict, out_list

def fill_with_none_and_cap(seq):
  if not len(seq): return []

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
  if addy == []: addy = [None] * len(model)
  i = 0
  while i < len(model):
    if model[i] == None and addy[i] != None:
      addy.insert(i, None)
    i+=1
  while len(addy) > len(model):
    addy.pop()

#
def pair_up(seqs0, seqs1):

  s0, l0, r0, extra_info0 = seqs0[:4]
  s1, l1, r1, extra_info1 = seqs1[:4]

  s0 = fill_with_none_and_cap(s0)
  s1 = fill_with_none_and_cap(s1)

  # print (s0)
  # print (s1)

  if s0 == []:
    # ???
    first0 = 1
    first1 = 1
  elif s1 == []:
    # ???
    first0 = 1
    first1 = 1
  elif extra_info0[0] == "BTN" and extra_info1[0] == "BTN":
    # simultaneos starters
    # find which element should be the first in the sequence
    init0, seqno0 = extra_info0[1], extra_info0[2]
    init1, seqno1 = extra_info1[1], extra_info1[2]
    if seqno0 > seqno1:
      # print ("1 initiator", init0, seqno0)
      first0 = seqno0 + 1
      first1 = init0
    else:
      # print ("0 initiator", init1, seqno1)
      first0 = init1
      first1 = seqno1

  elif extra_info0[0] == "BTN":
    # 0 is starter
    # s1 = [None] * (s1[0] - 1) + s1
    init1, seqno1 = extra_info1[0], extra_info1[1]
    first0 = 1
    first1 = seqno1
  elif extra_info1[0] == "BTN":
    # 1 is starter
    # s0 = [None] * (s0[0] - 1) + s0
    init0, seqno0 = extra_info0[0], extra_info0[1]
    first0 = seqno0 + 1
    first1 = 1
  else:
    print ("Not expected")
    raise Exception("BTN expected")

  # print first0, first1
  # at least one of them must start at 1
  while first0 > 1 and first1 > 1:
    first0 -= 1
    first1 -= 1
  # print first0, first1

  if s0 == [] and s1 == []:
    return (s0, l0, r0), (s1, l1, r1)

  if s0 == []:
    s0 = [None] * s1[-1]
  elif first0 > s0[0]:
    # trim sequence
    s0 = lfilter(lambda x: x is None or x>=first0, s0)
  elif first0 < s0[0]:
    # add trailling None
    s0 = [None] * (s0[0] - first0) + s0

  if s1 == []:
    s1 = [None] * s0[-1]
  elif first1 > s1[0]:
    # trim sequence
    s1 = lfilter(lambda x: x is None or x>=first1, s1)
  elif first1 < s1[0]:
    # add trailling None
    s1 = [None] * (s1[0] - first1) + s1

  # print (s0)
  # print (s1)
  # print()

  insert_nones(l0, s0)
  insert_nones(l1, s1)

  insert_nones(r0, s0)
  insert_nones(r1, s1)

  # print (zip (s0, s1)[-10:], len(s0), len(s1))
  # print (l0[-10:], len(l0))
  # print (l1[-10:], len(l1))

  # return lzip(s0, s1), zip(l0, l1), zip(r0, r1)
  return (s0, l0, r0), (s1, l1, r1)

def test_pair_up():

  # both nodes had button click before receiving a message
  n0 = [2, 4]
  n1 = [5]
  (n0, _, _), (n1, _, _) = pair_up((n0,n0,n0,("BTN",4,2)), (n1,n1,n1,("BTN",4,5)))
  print (lzip(n0, n1) == [ (None, None), (2, None), (None, None), (4, 5) ])

  n0 = [5]
  n1 = [2, 4]
  (n0, _, _), (n1, _, _) = pair_up((n0,n0,n0,("BTN",4,5)), (n1,n1,n1,("BTN",4,2)))
  print (lzip(n0, n1) == [ (None, None), (None, 2), (5, None) ])

  n0 = [2, 3]
  n1 = [100, 101]
  (n0, _, _), (n1, _, _) = pair_up((n0,n0,n0,("BTN",101,2)), (n1,n1,n1,("BTN",2,100)))
  print (lzip(n0, n1) == [ (None, None), (2, 100), (3, 101)])

  n0 = [1, 3]
  n1 = [3, 4, 6]
  (n0, _, _), (n1, _, _) = pair_up((n0,n0,n0,("BTN",4,1)), (n1,n1,n1,("BTN",1,3)))
  print (lzip(n0, n1) == [ (1, 3), (None, 4), (3, None) ])

  # only one node had button click
  n0 = [1, 3]
  n1 = [3, 4, 6]
  (n0, _, _), (n1, _, _) = pair_up((n0,n0,n0,("BTN",4,1)), (n1,n1,n1,(1,3,"BTN")))
  print (lzip(n0, n1) == [ (1, 3), (None, 4), (3, None) ])

  n0 = [3, 4, 6]
  n1 = [1, 3, 4]
  (n1, _, _), (n0, _, _) = pair_up((n0,n0,n0,(1,3,"BTN")), (n1,n1,n1,("BTN",4,1)))
  # print zip(n1, n0) == [ (3, 1), (4, None), (None, 3), (6, 4) ] # Old test
  print (lzip(n1, n0) == [ (4, 1), (None, None), (6, 3) ]) # test shifting to have always 0 as first sucessful deliver

  n0 = [8, 10, 12, 13]
  n1 = lrange(1,6)
  (n1, _, _), (n0, _, _) = pair_up((n0,n0,n0,(1,8)), (n1,n1,n1,("BTN",9,1)))
  print (lzip(n1, n0) == [ (None, 1), (10, 2), (None, 3), (12, 4), (13, 5) ])

  # unidirectional
  n0 = [4, 7]
  n1 = []
  (n1, _, _), (n0, _, _) = pair_up((n0,n0,n0,tuple()), (n1,n1,n1,tuple()))
  print (lzip(n1, n0) == [(None, None), (None, None), (None, None), (4, None), (None, None), (None, None), (7, None)])

def test_integrals():
  p1 = (0,0.0)
  p2 = (1,1.0)
  p3 = (0,1.0)
  p4 = (1,2.0)
  print (segment_diff_integral(p1, p2, p3, p4) == 1)
  p1 = (0,0.0)
  p2 = (1,1.0)
  p3 = (0,1.0)
  p4 = (1,10.0)
  print (segment_diff_integral(p1, p2, p3, p4) == 5)
  p1 = (3,0.0)
  p2 = (4,1.0)
  p3 = (3,1.0)
  p4 = (4,10.0)
  print (segment_diff_integral(p1, p2, p3, p4) == 5)
  p1 = (0,1.0)
  p2 = (1,0.0)
  p3 = (0,0.0)
  p4 = (1,1.0)
  print (segment_diff_integral(p1, p2, p3, p4) == 0.5)
  p1 = (3,2.0)
  p2 = (4,1.0)
  p3 = (3,1.0)
  p4 = (4,2.0)
  print (segment_diff_integral(p1, p2, p3, p4) == 0.5)
  p1 = (3,2.0)
  p2 = (4,1.0)
  p3 = (3,1.0)
  p4 = (4,11.0)
  print (int(100*segment_diff_integral(p1, p2, p3, p4)) == 459)
  p1 = (3,2.0)
  p2 = (5,1.0)
  p3 = (3,1.0)
  p4 = (5,2.0)
  print (segment_diff_integral(p1, p2, p3, p4) == 1)
  print (integreate_sequence_of_segments(range(0, 10, 2),range(1,10, 2))==1)
  print (integreate_sequence_of_segments(range(0, 10, 2),range(2,12, 2))==10)


def test_asimmetry_metric():
  mylen = 4
  metrics = []
  for i in (3,): #lrange(2**mylen):
    for j in (12,): #lrange(2**mylen):
      print ("i,j", i,j)
      s_test_0 = lrange(1, mylen+1)
      s_test_1 = lrange(1, mylen+1)
      for k in lrange(mylen):
        if i & 2**k:
          s_test_0[k] = None
        if j & 2**k:
          s_test_1[k] = None
      m = asimmetry_metric_kpe(s_test_0, s_test_1)
      metrics += [m]
      print(m, asimmetry_metric_simplified(s_test_0, s_test_1))
  # print (metrics)
  # return
  for i in sorted(set(metrics)):
    print ("%3.8f" % (i,), metrics.count(i))

  mylens = [8, 10, 12, 14]
  mylens = [7, 9, 11, 13, 15, 100,]# 500, 1000]
  # mylens = [3600]
  metrics = []
  for l in mylens:
    s_test_0 = [None]*(l//2) + lrange(l//2 + 1, l+1)
    s_test_1 = lrange(1, (l+1)//2 + 1) + [None]*(l//2)
    m = asimmetry_metric_kpe(s_test_0, s_test_1)
    metrics += [ m ]
    print(m, asimmetry_metric_simplified(s_test_0, s_test_1))

    # Random test
    # s_test_0 = []
    # s_test_1 = []
    # count = 1
    # while len(s_test_0) < l:
    #   if random.random() > 0.5:
    #     s_test_0 += [None]
    #   else:
    #     s_test_0 += [count]
    #
    #   if random.random() > 0.5:
    #     s_test_1 += [None]
    #   else:
    #     s_test_1 += [count]
    #
    #   count += 1
    # m = asimmetry_metric_kpe(s_test_0, s_test_1)
    # metrics += [ m ]

    # s_test_0 = [None]*((l+1)//2) + lrange((l+1)//2 + 1, l+1)
    # s_test_1 = lrange(1, (l)//2 + 1) + [None]*((l+1)//2)
    # print (s_test_0)
    # print (s_test_1)
  # print metrics
  for i in metrics:
    print (i)

def test_interpolate():
  test_dict = {}
  for i in lrange(5):
    test_dict[i] = i
  for i in lrange(15, 20):
    test_dict[i] = 1

  test_dict.pop(0)
  test_dict.pop(16)

  # for k,v in sorted(test_dict.items()):
  #   print (k,v)

  test_dict, test_list = interpolate(test_dict, interval=1, window=2, avg_method="Even")

  # for k,v in sorted(test_dict.items()):
  #   print (k,v)

  print (test_list == [1.5, 1.0, 2.0, 3.0, 4.0, 2.25, 2.25, 2.25, 2.25, 2.25,
                       2.25, 2.25, 2.25, 2.25, 2.25, 1.0, 1.75, 1.0, 1.0, 1.0])
################################################################################

def parse_twofiles(fname0, fname1, n_data_points=None, should_print = False):
  should_plot = False
  # fname0 = "node0.log"
  # fname1 = "node1.log"
  n0 = parsefile(fname0)
  n1 = parsefile(fname1)

  # print (n0[0],n0[-1])
  # print (n1[0],n1[-1])

  n0, n1 = pair_up(n0, n1)
  (s0, l0, r0), (s1, l1, r1) = n0, n1
  # print (lmap(len, (n0 + n1)), len(n0))
  # print (lzip(s0,s1)[:10], len(s0))

  if n_data_points is not None:
    s0 = s0[:min(len(s0),n_data_points)] + [None] * (n_data_points - len(s0))
    l0 = l0[:min(len(l0),n_data_points)] + [None] * (n_data_points - len(l0))
    r0 = r0[:min(len(r0),n_data_points)] + [None] * (n_data_points - len(r0))
    s1 = s1[:min(len(s1),n_data_points)] + [None] * (n_data_points - len(s1))
    l1 = l1[:min(len(l1),n_data_points)] + [None] * (n_data_points - len(l1))
    r1 = r1[:min(len(r1),n_data_points)] + [None] * (n_data_points - len(r1))

  ############
  # Delivery #
  ############
  pdr0 = 100.0*len(lfilter(lambda x: x!=None, s0))/len(s0)
  pdr1 = 100.0*len(lfilter(lambda x: x!=None, s1))/len(s1)
  if should_print: print ("Delivery rate 0:", pdr0)
  if should_print: print ("Delivery rate 1:", pdr1)

  loss_streak_avg0, loss_streak_stdev0 = average_loss_streak(s0)
  if should_print: print ("average_loss_streak(s0):", loss_streak_avg0, loss_streak_stdev0)
  loss_streak_avg1, loss_streak_stdev1 = average_loss_streak(s1)
  if should_print: print ("average_loss_streak(s1):", loss_streak_avg1, loss_streak_stdev1)
  delivery_streak_avg0, delivery_streak_stdev0 = average_delivery_streak(s0)
  if should_print: print ("average_delivery_streak(s0):", delivery_streak_avg0, delivery_streak_stdev0)
  delivery_streak_avg1, delivery_streak_stdev1 = average_delivery_streak(s1)
  if should_print: print ("average_delivery_streak(s1):", delivery_streak_avg1, delivery_streak_stdev1)
  asimmetry_metric = -1
  # asimmetry_metric = asimmetry_metric_kpe(s0, s1)
  if should_print: print ("asimmetry_metric based on expected instant probabilities:", asimmetry_metric)
  asimmetry_metric2 = asimmetry_metric_simplified(s0, s1)
  if should_print: print ("asimmetry_metric simplified:", asimmetry_metric2)

  if should_plot:
    plt.plot(lrange(len(s0)), lmap(lambda x: None if x == None else 1.1, s0), label="s0", marker='o', color = "blue")
    plt.plot(lrange(len(s0)), lmap(lambda x: 1.1 if x == None else None, s0), label="s0", marker='o', color = "red")
    plt.plot([ x+0.5 for x in lrange(len(s1))], lmap(lambda x: None if x == None else 0.9, s1), label="s1", marker='o', color = "green")
    plt.plot([ x+0.5 for x in lrange(len(s1))], lmap(lambda x: 0.9 if x == None else None, s1), label="s1", marker='o', color = "red")
    # plt.ylim( (0, 1.1) )
    plt.legend()
    plt.show()

  #######
  # LQI #
  #######
  l0_time_dict = {}
  for i, v in lzip(lrange(len(l0)), l0):
    if not v is None: l0_time_dict[i*2] = v

  l1_time_dict = {}
  for i, v in lzip(lrange(len(l1)), l1):
    if not v is None: l1_time_dict[1 + i*2] = v

  l0_interpol_dict, l0_interpol_list = interpolate(l0_time_dict, window=1, avg_method="Linear")
  l1_interpol_dict, l1_interpol_list = interpolate(l1_time_dict, window=1, avg_method="Linear")
  l0_interpol_list = l0_interpol_list[:min(len(l0_interpol_list), len(l1_interpol_list))]
  l1_interpol_list = l1_interpol_list[:min(len(l0_interpol_list), len(l1_interpol_list))]

  if len(l0_interpol_list) > 2 and len(l1_interpol_list) > 2:
    l_pearson, l_pearson_p = scipy.stats.pearsonr(l0_interpol_list, l1_interpol_list)
  else:
    l_pearson, l_pearson_p = None, None
  if should_print: print ("LQI pearson index:", l_pearson, l_pearson_p)

  l0_avg, l0_stdev = lmap(lambda x: (np.mean(x), np.std(x)), [lfilter(lambda x: x != None, l0)])[0]
  l1_avg, l1_stdev = lmap(lambda x: (np.mean(x), np.std(x)), [lfilter(lambda x: x != None, l1)])[0]
  if should_print: print ("LQI node 0 avg, stdev", l0_avg, l0_stdev)
  if should_print: print ("LQI node 1 avg, stdev", l1_avg, l1_stdev)

  # time x LQI plot
  if should_plot:
    plt.plot(lrange(0, 2*len(l0), 2), l0, label="l0", marker='o')
    plt.plot(lrange(len(l0_interpol_list)), l0_interpol_list, label="l0_inter", marker='+')
    plt.plot(lrange(1, 2*len(l1), 2), l1, label="l1", marker='o')
    plt.plot(lrange(len(l1_interpol_list)), l1_interpol_list, label="l1_inter", marker='+')
    plt.xlabel('t')
    plt.ylabel('LQI')
    plt.legend(loc='best')
    plt.show()
    plt.scatter(l0_interpol_list, l1_interpol_list, marker='o')
    plt.xlabel('l0')
    plt.ylabel('l1')
    plt.show()

  ########
  # RSSI #
  ########
  r0_time_dict = {}
  for i, v in lzip(lrange(len(r0)), r0):
    if not v is None: r0_time_dict[i*2] = v

  r1_time_dict = {}
  for i, v in lzip(lrange(len(r1)), r1):
    if not v is None: r1_time_dict[1 + i*2] = v

  r0_interpol_dict, r0_interpol_list = interpolate(r0_time_dict, window=1, avg_method="Linear")
  r1_interpol_dict, r1_interpol_list = interpolate(r1_time_dict, window=1, avg_method="Linear")
  r0_interpol_list = r0_interpol_list[:min(len(r0_interpol_list), len(r1_interpol_list))]
  r1_interpol_list = r1_interpol_list[:min(len(r0_interpol_list), len(r1_interpol_list))]

  if len(r0_interpol_list) > 2 and len(r1_interpol_list) > 2:
    r_pearson, r_pearson_p = scipy.stats.pearsonr(r0_interpol_list, r1_interpol_list)
  else:
    r_pearson, r_pearson_p = None, None
  if should_print: print ("RSSI pearson index:", r_pearson, r_pearson_p)

  r0_avg, r0_stdev = lmap(lambda x: (np.mean(x), np.std(x)), [lfilter(lambda x: x != None, r0)])[0]
  r1_avg, r1_stdev = lmap(lambda x: (np.mean(x), np.std(x)), [lfilter(lambda x: x != None, r1)])[0]
  if should_print: print ("RSSI node 0 avg, stdev", r0_avg, r0_stdev)
  if should_print: print ("RSSI node 1 avg, stdev", r1_avg, r1_stdev)

  # time x RSSI plot
  if should_plot:
    plt.plot(lrange(0, 2*len(r0), 2), r0, label="r0", marker='o')
    plt.plot(lrange(len(r0_interpol_list)), r0_interpol_list, label="r0_inter", marker='+')
    plt.plot(lrange(1, 2*len(r0), 2), r1, label="r1", marker='o')
    plt.plot(lrange(len(r1_interpol_list)), r1_interpol_list, label="r1_inter", marker='+')
    plt.xlabel('t')
    plt.ylabel('RSSI')
    plt.legend(loc="best")
    plt.show()
    plt.scatter(r0_interpol_list, r1_interpol_list, marker='o')
    plt.xlabel('r0')
    plt.ylabel('r1')
    plt.show()

  return (len(s0), asimmetry_metric, asimmetry_metric2, pdr0, pdr1, abs(pdr0-pdr1),
          delivery_streak_avg0, delivery_streak_avg1,
          delivery_streak_stdev0, delivery_streak_stdev1,
          loss_streak_avg0, loss_streak_avg1,
          loss_streak_stdev0, loss_streak_stdev1,
          l0_avg, l1_avg, l0_stdev, l1_stdev, l_pearson, l_pearson_p,
          r0_avg, r1_avg, r0_stdev, r1_stdev, r_pearson, r_pearson_p)




if __name__ == "__main__":
  # test_pair_up()
  # test_asimmetry_metric()
  # test_interpolate()
  # test_integrals()
  #
  # exit()

  if len(sys.argv) > 1:
    for f in sys.argv[1:]:
      try:
        print ("\nProcessing:", f + "node0.log", f + "node1.log")
        result = parse_twofiles(f + "node0.log", f + "node1.log", n_data_points=36000, should_print = True)
      except IOError as e:
        print ("Could not open file:", e, file=sys.stderr)
    exit()

  distance = ("moderate", "close", "far")
  power = ("1x1", "1x2", "3x3", "3x4", "7x7") # PA-LEVEL, value set at cc2420_set_txpower()
  relative_position = ("BxT", "TxB", "LxR")
  location = ("indoor", )#"outdoor")
  reps = lmap(str, range(1,2))
  factors = (distance, power, relative_position, location, reps)
  identifiers = ("D", "P", "R", "L", "i")
  identifiers_labels = {"D":"distance", "P":"power", "R":"relative position",
                        "L":"location", "i":"i"}

  def mix(l, it=0, res=[""]):
    if not len(l):
      return res
    newres = []
    for i in l[0]:
      for r in res:
        newres += [r+identifiers[it]+"-"+i+"_"]
    return mix(l[1:], it+1, newres)

  # filename format: D-${D}_P-${P}_R-${R}_L-${L}_i-${i}_node${node}.log
  files_prefixes = sorted(mix(factors))
  for i in files_prefixes:
    print ("\t",i)
  print ()

  out_file = open("_link_test_result.csv", 'w', encoding="ascii")

  headers = tuple(map(identifiers_labels.get, identifiers)) + (
          "# points", "Assimetry metric", "Assimetry metric simplified", "Avg delivery 0",
          "Avg delivery 1", "Avg delivery diff",
          "Avg delivery streak 0", "Avg delivery streak 1",
          "Stdev delivery streak 0", "Stdev delivery streak 1",
          "Avg loss streak 0", "Avg loss streak 1",
          "Stdev loss streak 0", "Stdev loss streak 1",
          "LQI avg 0", "LQI avg 1", "LQI stdev 0", "LQI stdev 1",
          "LQI pearson", "LQI pearson p", "RSSI avg 0", "RSSI avg 1",
          "RSSI stdev 0", "RSSI stdev 1", "RSSI pearson", "RSSI pearson p")

  print (";".join(headers))
  out_file.write(u";".join(headers) + "\n")

  for prefix in files_prefixes:
    try:
      print ("Processing:", prefix + "node0.log", prefix + "node1.log")
      result = parse_twofiles(prefix + "node0.log", prefix + "node1.log", n_data_points=36000)
      print (";".join(map(str,result)))
      out_file.write(prefix.replace("_",";") + u";".join(map(str,result)) + "\n")
    except IOError as e:
      print ("Could not open file:", e, file=sys.stderr)

  out_file.close()
