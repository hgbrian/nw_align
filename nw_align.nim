import tables
import sequtils

const
  NO = 0
  UP = 1
  LF = 2
  DG = 3
  ntlist = ["","A","C","G","T","-"]
  ntdict = {'A':1, 'C':2, 'G':3, 'T':4, '-':5}.toTable
  eps = 1e-8

# -------------------------------------------------------------
# Convert between dna sequence string and array for computation
#
proc unconvert_dna(align: seq): string {.inline.} =
  var
    n = len(align)-1
    align_str = ""
  while n >= 0:
    align_str.add(ntlist[align[n]])
    n -= 1
  return align_str


proc convert_dna(dnaseq: string): seq {.inline.} =
  let len_seq = len(dnaseq)
  var arr = repeat(0, len_seq)
  for ix in 0..<len_seq: arr[ix] = ntdict[dnaseq[ix]]
  return arr


proc sw_align(seqi: string, seqj: string, gap: float, match: float, mismatch: float):
              tuple[align_i: string, align_j: string] =
  # ------------------------------------------
  # Allocate seq arrays and turn chars to ints
  #

  let max_i = len(seqi)
  let max_j = len(seqj)

  var arri = convert_dna(seqi)
  var arrj = convert_dna(seqj)

  # ---------------------------------
  # Set up 2D arrays, score and point
  #
  var score = newSeqWith(max_i+1, newSeq[float](max_j+1))
  var point = newSeqWith(max_i+1, newSeq[int](max_j+1))


  for j in 1..max_j:
    point[0][j] = LF
    score[0][j] = 0

  for i in 1..max_i:
    point[i][0] = UP
    score[i][0] = 0

  score[0][0] = 0.0
  point[0][0] = NO


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

      # Smith-Waterman score never goes below 0
      dg_score = max(0, score[i-1][j-1] + matchscore)
      up_score = max(0, score[i-1][j] + gap)
      lf_score = max(0, score[i][j-1] + gap)

      if dg_score >= up_score and dg_score >= lf_score:
          score[i][j] = dg_score
          point[i][j] = DG
      elif lf_score >= dg_score and lf_score >= up_score:
          score[i][j] = lf_score
          point[i][j] = LF
      else:
          score[i][j] = up_score
          point[i][j] = UP


  # --------------------------------------
  # Find the max score in the entire array
  #
  var max_score = -1.0
  var max_score_pos = @[-1,-1]

  for i in 0..max_i:
    for j in 0..max_j:
      if score[i][j] > max_score:
        max_score = score[i][j]
        max_score_pos = @[i,j]

  echo $max_score, " ", $max_score_pos

  # ---------------------------------
  # Traverse
  #
  var
    align_i = repeat(0, max_i + max_j)
    align_j = repeat(0, max_i + max_j)
    nij = 0
    i = max_score_pos[0]
    j = max_score_pos[1]
    p = point[i][j]
    s = score[i][j]

  while s > eps:
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

    p = point[i][j]
    s = score[i][j]
    nij += 1


  # -------------------------
  # Output reversed sequence
  # Return a tuple of results
  #
  let ret: tuple[align_i: string, align_j: string] = (unconvert_dna(align_i), unconvert_dna(align_j))
  return ret


proc nw_align(seqi: string, seqj: string, gap: float, match: float, mismatch: float):
              tuple[align_i: string, align_j: string] =
  # ------------------------------------------
  # Allocate seq arrays and turn chars to ints
  #

  let max_i = len(seqi)
  let max_j = len(seqj)

  var arri = convert_dna(seqi)
  var arrj = convert_dna(seqj)

  # ---------------------------------
  # Set up 2D arrays, score and point
  #
  var score = newSeqWith(max_i+1, newSeq[float](max_j+1))
  var point = newSeqWith(max_i+1, newSeq[int](max_j+1))

  for j in 1..max_j:
    point[0][j] = LF
    score[0][j] = gap * float(j)

  for i in 1..max_i:
    point[i][0] = UP
    score[i][0] = gap * float(i)

  score[0][0] = 0.0
  point[0][0] = NO
  echo "score[0][0] ", $score[0][0]

  # ---------------------------------
  # Scores
  #
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

      dg_score = score[i-1][j-1] + matchscore
      up_score = score[i-1][j] + gap
      lf_score = score[i][j-1] + gap
  
      if dg_score >= up_score and dg_score >= lf_score:
          score[i][j] = dg_score
          point[i][j] = DG
      elif lf_score >= dg_score and lf_score >= up_score:
          score[i][j] = lf_score
          point[i][j] = LF
      else:
          score[i][j] = up_score
          point[i][j] = UP

  # ---------------------------------
  # Traverse
  #
  var
    align_i = repeat(0, max_i + max_j)
    align_j = repeat(0, max_i + max_j)
    nij = 0
    i = max_i
    j = max_j
    p = point[i][j]
    s = score[i][j]

  while i > 0 or j > 0:
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

    p = point[i][j]
    s = score[i][j]
    nij += 1


  # -------------------------
  # Output reversed sequence
  # Return a tuple of results
  #
  let ret: tuple[align_i: string, align_j: string] = (unconvert_dna(align_i), unconvert_dna(align_j))
  return ret

  
# ----------------------------------------------------------------------------------------
# Get sequences seqi, seqj as input
# params
#
let
  gap = -1.0
  match = float(1.0)
  mismatch = -1.0


when not defined(js): # defined(c) not working?
  import os
  let
    kind = paramStr(1)
    seqi = paramStr(2)
    seqj = paramStr(3)

  if kind == "local":
    var alignment = sw_align($seqi, $seqj, gap, match, mismatch)
    echo alignment.align_i & "\n" & alignment.align_j
  elif kind == "global":
    var alignment = nw_align($seqi, $seqj, gap, match, mismatch)
    echo alignment.align_i & "\n" & alignment.align_j
  else:
    echo "alignment type must be local or global"


when defined(js):
  # add console.log(nw_align_js(process.argv[2], process.argv[3])); to nw_align.js
  # as the last line for command-line javascript version
  proc nw_align_js(seqi: cstring, seqj: cstring): cstring {.exportc.} =
    var alignment = nw_align($seqi, $seqj, gap, match, mismatch)
    result = alignment.align_i & "\n" & alignment.align_j

  proc sw_align_js(seqi: cstring, seqj: cstring): cstring {.exportc.} =
    var alignment = sw_align($seqi, $seqj, gap, match, mismatch)
    result = alignment.align_i & "\n" & alignment.align_j
