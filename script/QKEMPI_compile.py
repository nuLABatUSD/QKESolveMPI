import os
import sys

prev_res = []
for f in os.listdir("results/"):
    poss = f.split('-')[0]
    if poss.isnumeric():
        prev_res.append(int(poss))

if len(prev_res) == 0:
    next = 1
else:
    next = max(prev_res)+1

filename = sys.argv[1].split('/')[-1]
filename = filename.split('.')[0]

if next < 10:
    print("0{}-{}".format(next, filename))
else:
    print("{}-{}".format(next, filename))
