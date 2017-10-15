import subprocess
from timeit import default_timer as timer
from random import choice

out = {}

cmds = {
        "orig":      ["python2", "nw_align_original.py"],
        "orig3":     ["python", "nw_align_original.py"],
        "numpy":     ["python", "nw_align.py", "numpy"],
        "numba":     ["python", "nw_align.py", "numba"],
        "torch":     ["python", "nw_align.py", "torch"],
        "torchcuda": ["python", "nw_align.py", "torchcuda"],
        "cupy":      ["python", "nw_align.py", "cupy"],
        "nimc":      ["./nw_align"],
        "nimjs":     ["node", "nimcache/nw_align.js"],
        "js":        ["node", "nw_align.js"]
       }

Ns = [500,1000,1500,2000,2500,3000,5000]
maxlen = {"cupy":1000, "torchcuda":5000}

out = {}
out["out"] = {k:{} for k in cmds}
out["time"] = {k:{} for k in cmds}

for N in Ns:
  seq1 = ''.join(choice("ACGT") for _ in range(N))
  seq2 = ''.join(choice("ACGT") for _ in range(N))

  for k, cmd in cmds.items():
    if k in maxlen and N > maxlen[k]: continue

    start = timer()
    out["out"][k][N] = subprocess.check_output(cmd + [seq1, seq2]).decode()
    end = timer()
    out["time"][k][N] = end - start

for N in Ns:
  outs = {k: out["out"][k][N] for k in cmds}
  assert len(set(outs.values())) == 1, "not the same output from every method {}".format(set(outs.values()))
  # DEBUG print(outs)

print(out["time"])
