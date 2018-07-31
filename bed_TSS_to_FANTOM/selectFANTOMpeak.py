import sys

d = {}
infile = open(sys.argv[1], 'r')
for line in infile:

     line = line.strip().split("\t")
     ID = line[3]

     if ID in d:
          d[ID].append(line)
     else:
          d[ID] = [line]

for ID in d:

     d[ID].sort(key=lambda x: int(x[4]), reverse=True)
     print("\t".join(d[ID][0]))






