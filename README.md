Needleman Wunsch
----------------
A basic no-frills implementation.

Testing Nim vs Python vs Javascript

    wget https://bitbucket.org/brentp/biostuff/raw/282b504ac9020fe1449e23f800b20b5bd7d12061/nwalign/pairwise.py -O nw_align_original.py
    sed -i s/print\ global_align\(sys.argv\\[1\\],\ sys.argv\\[2\\]\)/print\(global_align\(sys.argv\[1\],\ sys.argv\[2\]\)/ nw_align_original.py
    sed -i s/print\ read_matrix/#print\ read_matrix/ nw_align_original.py
    nim c -d:release nw_align.nim
    nim js -d:release nw_align.nim
    echo "console.log(nw_align_js(process.argv[2], process.argv[3]));" >>nimcache/nw_align.js
    python test_nw_align.py
