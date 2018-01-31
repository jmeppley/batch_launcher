#!/usr/bin/env bats

@test "Testing HMM db split" {
    run ./batch_launcher.py -X local -N 5 -- hmmsearch --tblout test_data/.tst.genes.v.pfam.tbl --cpu 2 test_data/pfams.100.hmm test_data/genes.faa
    [ "$status" = 0 ]
    run diff <(grep -v "^#" test_data/.tst.genes.v.pfam.tbl) <(grep -v "^#" test_data/genes.v.pfam.tbl)
    [ "$status" = 0 ]
}

@test "Testing HMM db split with keep" {
    run ./batch_launcher.py -X local -N 5 -K -- hmmsearch --tblout test_data/.tst.genes.v.pfam.tbl --cpu 2 test_data/pfams.100.hmm test_data/genes.faa
    [ "$status" = 0 ]
    run diff <(grep -v "^#" test_data/.tst.genes.v.pfam.tbl) <(grep -v "^#" test_data/genes.v.pfam.tbl)
    [ "$status" = 0 ]
}

