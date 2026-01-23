import os
import sys

prev_res = []
for f in os.listdir("results/"):
    if f[:2].isnumeric():
        prev_res.append(int(f[:2]))

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
