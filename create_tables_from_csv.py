from collections import defaultdict

def whiteToBlack(mix=0):
  mix = min(1, mix)
  mix = max(0, mix)
  mix = 255 - int(255*mix)
  return ("%02X" % (mix,)) * 3

f = file("_link_test_result.csv")

f.readline()

values = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: None))))

powers = set()
positions = set()
distances = set()

for line in f:
  (distance, power, position, location, iteration, points, assim_metric,
  assim_metric_simplified, avg_delivery0, avg_delivery1, avg_delivery_diff,
  avg_delivery_streak0, avg_delivery_streak1, _, _,
  avg_loss_streak0, avg_loss_streak1, _, _,
  LQI_avg0, LQI_avg1, _, _, LQI_pearson, LQI_pearson_p,
  RSSI_avg0, RSSI_avg1,  _, _, RSSI_pearson, RSSI_pearson_p
  ) = line.strip().split(";")

  position = position[1+position.find("-"):]
  distance = distance[1+distance.find("-"):]
  distance = str.upper(distance[0])

  values[power][position][distance] = [assim_metric, assim_metric_simplified]

  powers.add(power)
  positions.add(position)
  distances.add(distance)


powers = sorted(powers)
positions = list(positions)
# distances = sorted(distances)

for power in powers:
  fout = file("table_"+power+".tex", 'w')
  fout.write(r"""
\begin{table}[htb]
\centering
\begin{tabular}{cccccc}
\cline{3-5}
""")
  print power

  for position in positions:
    if position == positions[-1]:
      fout.write(r"""\multirow{-3}{*}{Positioning} &""")
    else:
      fout.write(r"""                              &""")
    fout.write(r"""\multicolumn{1}{c|}{""")
    fout.write(position)
    fout.write(r"""} & """)
    print position,
    for distance in distances:
      m = values[power][position][distance][1]
      if m != None:
        m = float(m)
        # print power, position, distance
        fout.write(r"""\multicolumn{1}{c|}{\cellcolor[HTML]{""")
        fout.write(whiteToBlack(m))
        fout.write(r"""}{\color[HTML]""")
        if m > 0.4:
          fout.write(r"{FFFFFF} \textbf{")
        else:
          fout.write("{000000}          ")
        fout.write("%.2f" % (m,))
        if m > 0.4: fout.write("}")
        fout.write(r""" }} & """)
      else:
        fout.write(r"""\multicolumn{1}{c|}{\cellcolor[HTML]{FFFFFF}{\color[HTML]{000000} - }} & """)
      print values[power][position][distance][1],
    if position == positions[-1]:
      fout.write(r"""\multirow{-3}{*}{\includegraphics[height=40pt]{rect.png}} \\ \cline{3-5}""")
    else:
      fout.write(r"""                   \\ \cline{3-5}""")
    fout.write("\n")
    print

  fout.write(r""" & & """)
  print "\t",
  for distance in distances:
    fout.write(distance)
    fout.write(r""" & """)
    print distance,
  fout.write(r"""\\
 & & \multicolumn{3}{c}{Distances} &
\end{tabular}
\caption{Asymmetry metric values -- """)
  fout.write(power.replace("P-", "Power combination "))
  fout.write(r"""}
\end{table}
""")
  print
  print
  fout.close()


  # print whiteToBlack(float(assim_metric_simplified)), assim_metric_simplified
