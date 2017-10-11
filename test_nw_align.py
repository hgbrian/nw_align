import subprocess
import time
from random import choice

out = {}

cmds = {"orig": ["python2", "nw_align_original.py"],
        "numpy": ["python", "nw_align.py", "numpy"],
        "numba": ["python", "nw_align.py", "numba"],
        #"torchcuda": ["python", "pairwise.py", "torchcuda"],
        "nim": ["./nw_align"],
        #"nimjs": ["nodejs", "nimcache/nw_align.js"]
        }

for N in range(500,5500,500):
  out[N] = {"out":{}, "time":{}}

  seq1 = ''.join(choice("ACGT") for _ in range(N))
  seq2 = ''.join(choice("ACGT") for _ in range(N))

  for k, cmd in cmds.items():
    t = time.time()
    out[N]["out"][k] = subprocess.check_output(cmd + [seq1, seq2])
    out[N]["time"][k] = time.time() - t


for N in out:
  print(N)
  for k,v in out[N]["out"].items():
    print("  {} {:.3g}".format(k, out[N]["time"][k]))
    #print("  {}".format(v.decode()))
