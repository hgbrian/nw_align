import tables
import sequtils

const
  NO = 0
  UP = 1
  LF = 2
  DG = 3
  ntlist = ["","A","C","G","T","-"]
  ntdict = {'A':1, 'C':2, 'G':3, 'T':4, '-':5}.toTable


proc align(seqi: string, seqj: string, gap: float, match: float, mismatch: float): 
           tuple[align_i: string, align_j: string] =
    # ------------------------------------------
    # Allocate seq arrays and turn chars to ints
    #

    let 
      max_i = len(seqi)
      max_j = len(seqj)

    var arri = repeat(0, max_i)
    var arrj = repeat(0, max_j)

    for ix in 0..max_i-1: arri[ix] = ntdict[seqi[ix]]
    for ix in 0..max_j-1: arrj[ix] = ntdict[seqj[ix]]


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

    var
      align_i = repeat(0, max_i + max_j)
      align_j = repeat(0, max_i + max_j)
      nij = 0
      p = 0
      s = 0.0

    var
      i = max_i
      j = max_j

    while i > 0 or j > 0:
      p = point[i][j]
      s = score[i][j]

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
    proc unconvert_dna(align: seq): string =
      var 
        n = len(align)-1
        align_str = ""
      while n >= 0:
        align_str.add(ntlist[align[n]])
        n -= 1
      return align_str
  

    # Return a tuple of results
    let ret: tuple[align_i: string, align_j: string] = (unconvert_dna(align_i), unconvert_dna(align_j))
    return ret

  
# ----------------------------------------------------------------------------------------
# Get sequences seqi, seqj as input
#
# params
let
  gap = -1.0
  match = float(1.0)
  mismatch = -1.0
  

when not defined(js): # defined(c) not working?
  import os
  let 
    seqi = paramStr(1)
    seqj = paramStr(2)

  var alignment = align($seqi, $seqj, gap, match, mismatch)
  echo alignment.align_i & "\n" & alignment.align_j


when defined(js):
  proc nw_align_js(seqi: cstring, seqj: cstring): cstring {.exportc.} =
    var alignment = align($seqi, $seqj, gap, match, mismatch)
    result = alignment.align_i & "\n" & alignment.align_j

