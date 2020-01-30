import operator

def check():
  seqno_sequence = filter(None, seqno_sequence)
  print seqno_sequence
  i = 1
  while i < len(seqno_sequence):
    if seqno_sequence[i] > seqno_sequence[i - 1] + 1:
      missed_packets = seqno_sequence[i] - seqno_sequence[i - 1] - 1
      print "miss", seqno_sequence[i-1], missed_packets
    i+= 1

def average_distance_to_nearest_paired_loss(s0, s1, swap):
  distances = []
  max_i = min(len(s0),len(s1))
  for i in range(max_i):
    # print "i", i
    if s0[i] == None:
      # print "None found"
      if s1[i] == None:
        distances += [0]
        continue
      distance = 1
      should_break = False
      while i + distance < max_i or i - distance >= 0:
        if i - distance >= 0 and s1[i - distance] == None:
          # print "before i", i, distance
          if not swap: distance -= 1
          should_break = True
        if i + distance < max_i and s1[i+distance] == None:
          # print "after i", i, distance
          if swap: distance -= 1
          should_break = True
        if should_break: break
        distance += 1
      if distance < max_i:
        distances += [distance]
      # else:
      #   print "out of bounds"

  return np.mean(distances), np.std(distances), distances

def test_average_distance_to_nearest_paired_loss():
  print ("test average_distance_to_nearest_paired_loss 1")
  s_test_0 = range(1, 10)
  s_test_0[5] = None
  s_test_1 = range(1, 10)
  s_test_1[2] = None
  s_test_1[6] = None
  d, _, _ = average_distance_to_nearest_paired_loss( s_test_0, s_test_1, False)
  print (d == 1)

  s_test_0 = range(1, 10)
  s_test_0[5] = None
  s_test_1 = range(1, 10)
  s_test_1[3] = None
  s_test_1[7] = None
  d, _, _ = average_distance_to_nearest_paired_loss( s_test_0, s_test_1, False)
  print (d == 1)

  print ("test average_distance_to_nearest_paired_loss 2")
  s_test_0 = range(1, 10)
  s_test_0[3] = None
  s_test_0[7] = None
  s_test_1 = range(1, 10)
  s_test_1[5] = None
  d, _, _ = average_distance_to_nearest_paired_loss( s_test_1, s_test_0, True)
  print (d == 1)

  s_test_0 = range(1, 10)
  s_test_0[4] = None
  s_test_0[8] = None
  s_test_1 = range(1, 10)
  s_test_1[5] = None
  d, _, _ = average_distance_to_nearest_paired_loss( s_test_1, s_test_0, True)
  print (d == 1)

  # for i in s_test_0:
  #   print "X" if i == None else i, "  ",
  # print
  # for i in s_test_1:
  #   print "  ", "X" if i == None else i,
  # print

def calculate_distances(s0, s1, swap):
  fail_distances = []
  ok_distances = []
  max_i = min(len(s0),len(s1))
  for i in range(max_i):
    # print "i", i
    if s0[i] == None:
      op = operator.eq
      distances = fail_distances
    else:
      op = operator.ne
      distances = ok_distances

    if op(s1[i], None):
      distances += [0]
      continue
    distance = 1
    should_break = False
    while i + distance < max_i or i - distance >= 0:
      if i - distance >= 0 and op(s1[i - distance], None):
        # print "before i", i, distance
        if not swap: distance -= 1
        should_break = True
      if i + distance < max_i and op(s1[i+distance], None):
        # print "after i", i, distance
        if swap: distance -= 1
        should_break = True
      if should_break: break
      distance += 1
    # while distance < max_i:
    #   if op(s1[i - distance], None):
    #     # print "before i", i, distance
    #     if not swap: distance -= 1
    #     should_break = True
    #   if op(s1[ (i+distance) % max_i], None):
    #     # print "after i", i, distance
    #     if swap: distance -= 1
    #     should_break = True
    #   if should_break: break
    #   distance += 1
    if not should_break:
      # print "not should_break", distance, i, max_i
      # distance = max_i - i + 1.0/max_i
      # distance = 1 + max_i/2.0 + min(i, max_i-i-1)
      distance = ((max_i-1)/2.0)
    else:
      distance = min(distance, max_i)
    distances += [distance]

  # print ok_distances
  # print fail_distances

  return ok_distances, fail_distances

def asimmetry_metric(s0, s1):
  should_print = False
  lists = []
  lists += calculate_distances(s0, s1, False)

  if should_print:
    i_fail = i_ok = 0
    for i in s0:
      if i == None:
        print "%5.1f" % (lists[1][i_fail],),
        i_fail += 1
      else:
        print "%5.1f" % (lists[0][i_ok],),
        i_ok += 1
    print
    print "   ",
    for i in s0:
      print "X" if i == None else i, "   ",
    print

  lists += calculate_distances(s1, s0, True)

  if should_print:
    i_fail = i_ok = 0
    print " "*6,
    for i in s1:
      print "X" if i == None else i, "   ",
    print
    print "  ",
    for i in s1:
      if i == None:
        print "%5.1f" % (lists[3][i_fail],),
        i_fail += 1
      else:
        print "%5.1f" % (lists[2][i_ok],),
        i_ok += 1
    print

    print lists

  lists = filter(lambda x: len(x), lists)
  means = map(lambda x: sum(x), lists)
  lens = map(lambda x: len(x), lists)

  m = 1.0*reduce(operator.add, means)/reduce(operator.add, lens)
  m /= ((min(len(s0),len(s1))-1)/2.0)
  if should_print:
    print m
    print
  return m

def kpe2(t, k_plus_vector, k_minus_vector):
  k_plus_t = 0
  k_minus_t = 0
  threshold = 1e-8

  t_int = min(len(k_plus_vector)-1, max(0, int(round(t))))

  for m, s in k_plus_vector[t_int:]:
    v = scipy.stats.norm(m, s).pdf(t)
    k_plus_t += v
    if v < threshold: break

  for m, s in reversed(k_plus_vector[:t_int]):
    v = scipy.stats.norm(m, s).pdf(t)
    k_plus_t += v
    if v < threshold: break

  for m, s in k_minus_vector[t_int:]:
    v = scipy.stats.norm(m, s).pdf(t)
    k_minus_t += v
    if v < threshold: break

  for m, s in reversed(k_minus_vector[:t_int]):
    v = scipy.stats.norm(m, s).pdf(t)
    k_minus_t += v
    if v < threshold: break

  kpe_t = k_plus_t/(k_plus_t + k_minus_t)
  return kpe_t

def kpe(t, k_plus_vector, k_minus_vector):
  k_plus_t = 0
  k_minus_t = 0
  for (m, s) in k_plus_vector:
    k_plus_t += scipy.stats.norm(m, s).pdf(t)
  for (m, s) in k_minus_vector:
    k_minus_t += scipy.stats.norm(m, s).pdf(t)
  kpe_t = k_plus_t/(k_plus_t + k_minus_t)
  # print "kpe_(%f)=" % (t,), kpe_t
  return kpe_t

def kpe_dict(t, k_plus_vector, k_minus_vector):
  k_plus_t = 0
  k_minus_t = 0

  for (m, s) in k_plus_vector:
    if (m,s,t) in kpe_dict_values:
      v = kpe_dict_values[(m, s, t)]
    else:
      v = scipy.stats.norm(m, s).pdf(t)
      kpe_dict_values[(m, s, t)] = v
    k_plus_t += v

  for (m, s) in k_minus_vector:
    if (m,s,t) in kpe_dict_values:
      v = kpe_dict_values[(m, s, t)]
    else:
      v = scipy.stats.norm(m, s).pdf(t)
      kpe_dict_values[(m, s, t)] = v
    k_minus_t += v

  kpe_t = k_plus_t/(k_plus_t + k_minus_t)
  return kpe_t


# m = scipy.integrate.quad(lambda t: abs(kpe_dict(t, s0_k_plus_gaussians, s0_k_minus_gaussians) - kpe_dict(t, s1_k_plus_gaussians, s1_k_minus_gaussians)), -1, max_i*2,  epsabs=max_err, epsrel=max_err)[0]/(max_i*2+2)
# if should_print:
#   print "m kpe_dict", m


# m = scipy.integrate.quad(lambda t: abs(kpe2(t, s0_k_plus_gaussians, s0_k_minus_gaussians) - kpe2(t, s1_k_plus_gaussians, s1_k_minus_gaussians)), -1, max_i*2,  epsabs=max_err, epsrel=max_err)[0]/(max_i*2+2)
# if should_print:
#   print "m kpe2", m

# m = scipy.integrate.quad(lambda t: abs(kpe(t, s0_k_plus_gaussians, s0_k_minus_gaussians) - kpe(t, s1_k_plus_gaussians, s1_k_minus_gaussians)), -1, max_i*2,  epsabs=max_err, epsrel=max_err)[0]/(max_i*2+2)
# if should_print:
#   print 0, 2*max_i-1, m

def unpair(s):
  s0, s1 = [], []
  for i in s:
    s0 += [i[0]]
    s1 += [i[1]]
  return s0, s1



# print average_distance_to_nearest_paired_loss(s0,s1,0)
# print average_distance_to_nearest_paired_loss(s1,s0,1)
# print "asimmetry_metric based on distances",
# print asimmetry_metric(s0,s1)

# LQI x LQI plot
# fig = plt.figure()
lqi_pairs = lzip(l0, l1)
plt.scatter(*zip(*lqi_pairs))
plt.xlabel('LQI node 0')
plt.ylabel('LQI node 1')
plt.show()

# LQI pearson coef
lqi_pairs = lfilter(lambda x: None not in x, lqi_pairs)
lqi_coef = np.corrcoef([x[0] for x in lqi_pairs], [x[1] for x in lqi_pairs])[0][1]
print ("pearson coef LQI", lqi_coef)

# RSSI x RSSI plot
# fig = plt.figure()
rssi_pairs = lzip(r0, r1)
plt.scatter(*lzip(*rssi_pairs))
plt.xlabel('RSSI node 0')
plt.ylabel('RSSI node 1')
plt.show()

# RSSI pearson coef
rssi_pairs = lfilter(lambda x: None not in x, rssi_pairs)
rssi_coef = np.corrcoef([x[0] for x in rssi_pairs], [x[1] for x in rssi_pairs])[0][1]
print ("pearson coef RSSI", rssi_coef)

  # files_prefixes = []

  # for i in lrange(len(factors)):
  #   print (factors[i])
  #   for k in lrange(min(1,i),len(factors[i])):
  #     mystr = ""
  #     for j in lrange(len(factors)):
  #       print (i,k,j)
  #       mystr += identifiers[j] + "-"
  #       if i != j:
  #         mystr += factors[j][0]
  #       else:
  #         mystr += factors[i][k]
  #       mystr += "_"
  #     print (mystr)
  #     for rep in reps:
  #       files_prefixes += [mystr + "i-" + rep + "_"]
  #   print("============")
