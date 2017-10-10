import os
import tables
import arraymancer
import sequtils

const
  NO = 0
  UP = 1
  LF = 2
  DG = 3
  
let
  gap = -1.0
  match = 1.0
  mismatch = -1.0
  
# ---------------------------------
# Get sequences seqi, seqj as input
#

let 
  seqi = paramStr(1)
  seqj = paramStr(2)

#let 
#  seqi = "ACGTGATG"
#  seqj = "GTGATGATG"

echo "seqs", seqi, " ", seqj

# ------------------------------
# nucleotide char to int
# "A" is a string, 'A' is a char
#

const ntlist = ["","A","C","G","T","-"]
const ntdict = {'A':1, 'C':2, 'G':3, 'T':4, '-':5}.toTable

echo ntdict['A'], " ", ntdict.len
assert ntdict.hasKey('C')
for key, value in ntdict:
  echo key, " " ,value


# -------------------------------------------
# Read in params
#


# ------------------------------------------
# Allocate seq arrays and turn chars to ints
#

let 
  max_i = len(seqi)
  max_j = len(seqj)

var arri = zeros([max_i], int)
var arrj = zeros([max_j], int)

echo "??? ", seqi[0], "A"
echo "??? ", ntdict['A']
echo "??? ", ntdict[seqi[0]]

for ix in 0..max_i-1: arri[ix] = ntdict[seqi[ix]]
for ix in 0..max_j-1: arrj[ix] = ntdict[seqj[ix]]

# ---------------------------------
# Set up 2D arrays, score and point
#
var score = zeros([max_i+1, max_j+1], float)
var point = zeros([max_i+1, max_j+1], int)

point[0,0] = NO
point[0,1.._] = LF
point[1.._,0] = UP

#let
#  nums = @[1, 2, 3, 4]
#  strings = nums.mapIt(string, $(4 * it))

#var a = toSeq(0..max_j)
#var b = mapIt(a, float)
#echo b

score[_, 0] = (gap * (toSeq(0..max_i).toTensor().astype(float)))
score[0, _] = (gap * (toSeq(0..max_j).toTensor().astype(float))).reshape(1, max_j+1)


var
  ci = 0
  cj = 0
  matchscore = 0.0
  dg_score = 0.0
  up_score = 0.0
  lf_score = 0.0

for i in 1..max_i:
  ci = arri[i-1]
  for j in 1..max_j:
    cj = arrj[j-1]

    if ci == cj:
      matchscore = match
    else:
      matchscore = mismatch

    dg_score = score[i-1, j-1] + matchscore
    up_score = score[i-1, j] + gap
    lf_score = score[i, j-1] + gap
    
    
    if dg_score >= up_score and dg_score >= lf_score:
        score[i, j] = dg_score
        point[i, j] = DG
    elif up_score >= dg_score and up_score >= lf_score:
        score[i, j] = up_score
        point[i, j] = UP
    else:
        score[i, j] = lf_score
        point[i, j] = LF

var
  align_i = zeros([max_i + max_j], int)
  align_j = zeros([max_i + max_j], int)
  nij = 0
  p = 0
  s = 0.0

var
  i = max_i
  j = max_j

while i > 0 or j > 0:
  p = point[i, j]
  s = score[i, j]

  if p == DG:
      align_j[nij] = arrj[j-1]
      align_i[nij] = arri[i-1]
      i -= 1
      j -= 1
  elif p == LF:
      align_j[nij] = arrj[j-1]
      align_i[nij] = ntdict['-']
      j -= 1
  elif p == UP:
      align_j[nij] = ntdict['-']
      align_i[nij] = arri[i-1]
      i -= 1
  else:
      echo "ERR", p
      break

  nij += 1


# ------------------------
# Output reversed sequence
#
var
  n = 0
  align_j_str = ""
  align_i_str = ""

n = align_j.shape[0]-1
while n >= 0:
  align_j_str.add(ntlist[align_j[n]])
  n -= 1
  
n = align_i.shape[0]-1
while n >= 0:
  align_i_str.add(ntlist[align_i[n]])
  n -= 1

echo align_j_str
echo align_i_str


