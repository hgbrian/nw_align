Needleman Wunsch
----------------
A basic no-frills implementation.

Testing Nim vs Python vs Javascript

    nim c -d:release nw_align.nim
    nim js -d:release nw_align.nim
    echo "console.log(nw_align_js(process.argv[2], process.argv[3]));" >>nimcache/nw_align.js
    python test_nw_align.py
