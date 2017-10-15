const
  UP = 0
  LF = 1
  DG = 2
  NO = 3;


// from string to uint8 list/array and back again
const ntlist = ["","A","C","G","T","-"]
const ntdict = {"":0, "A":1, "C":2, "G":3, "T":4, "-":5}
const eps = 0.0000001;

function _convert_dna(dnastr) {
    let dnalist = dnastr.split('');
    var dnai = [];
    for (var i=0; i<dnalist.length; i++) {
      dnai.push(ntdict[dnalist[i]]);
    }
    return dnai;
}

function _unconvert_dna(dnai) {
  var dnalist = [];
  for (var i=0; i<dnai.length; i++) {
    dnalist.push(ntlist[dnai[i]]);
  }
  return dnalist.join('');
}


function repeat(dimensions, val) {
    var array = [];

    for (var i = 0; i < dimensions[0]; ++i) {
        array.push(dimensions.length == 1 ? val : repeat(dimensions.slice(1), val));
    }

    return array;
}


function global_align(seqi, seqj, gap, match, mismatch) {
  if (gap===undefined) gap = -1;
  if (match===undefined) match = 1;
  if (mismatch===undefined) mismatch = -1;

  let arrj = _convert_dna(seqj)
  let arri = _convert_dna(seqi)

  let max_j = arrj.length;
  let max_i = arri.length;

  var score = repeat([max_i+1, max_j+1], 0);
  var point = repeat([max_i+1, max_j+1], 0);

  score[0][0] = 0.0;

  for (var j=1; j<max_j+1; j++) {
    point[0][j] = LF;
    score[0][j] = gap * j;
  }

  for (var i=1; i<max_i+1; i++) {
    point[i][0] = UP;
    score[i][0] = gap * i;
  }

  point[0][0] = NO;

  var
    ci = 0,
    cj = 0;

  for (var i=1; i<max_i+1; i++) {
    ci = arri[i-1];
    for (var j=1; j<max_j+1; j++) {
      cj = arrj[j-1];

      if (ci == cj) {
        matchscore = match;
      }
      else {
        matchscore = mismatch;
      }

      dg_score = score[i-1][j-1] + matchscore;
      up_score = score[i-1][j] + gap;
      lf_score = score[i][j-1] + gap;

      if (dg_score >= up_score-eps && dg_score >= lf_score-eps) {
        score[i][j] = dg_score;
        point[i][j] = DG;
      }
      else if (lf_score >= dg_score-eps && lf_score >= up_score-eps) {
        score[i][j] = lf_score;
        point[i][j] = LF;
      }
      else {
        score[i][j] = up_score;
        point[i][j] = UP;
      }
    }
  }

  var
    align_i = repeat(max_i + max_j, 0),
    align_j = repeat(max_i + max_j, 0),
    i = max_i,
    j = max_j,
    nij = 0,
    p = 0,
    s = 0.0;

  while (i > 0 || j > 0) {
    p = point[i][j];
    s = score[i][j];

    if (p == DG) {
      align_j[nij] = arrj[j-1];
      align_i[nij] = arri[i-1];
      i -= 1;
      j -= 1;
    }
    else if (p == LF) {
      align_j[nij] = arrj[j-1];
      align_i[nij] = ntdict['-'];
      j -= 1;
    }
    else if (p == UP) {
      align_j[nij] = ntdict['-'];
      align_i[nij] = arri[i-1];
      i -= 1;
    }
    else {
      console.log("ERR "+p);
      break
    }
    nij += 1;
  }

  return {"align_i":_unconvert_dna(align_i), "align_j":_unconvert_dna(align_j)}

}

function reverse_str(astr) {
  return astr.split("").reverse().join("");
}

function nw_align(seqi, seqj) {
  const x = 1;
  //console.log(seqi+" "+seqj);
  var dnai = _convert_dna(seqi);
  var dnaj = _convert_dna(seqj);

  var al = global_align(seqi, seqj);
  return reverse_str(al.align_i)+"\n"+reverse_str(al.align_j);
}

// node nw_align.js ACGT TGCA
var alstr = nw_align(process.argv[2], process.argv[3]);
console.log(alstr);
