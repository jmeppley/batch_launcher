#!/usr/bin/env bats

@test "Testing missing chunk flag" {
    run ./batch_launcher.py -X local -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 2 ]
}

@test "Testing missing wait flag" {
    run ./batch_launcher.py -X sge -C 10 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 2 ]
}

@test "Testing missing input designation" {
    run ./batch_launcher.py -C 10 -X local -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 1 ]
}

@test "Testing wrong input designation" {
    run ./batch_launcher.py -C 10 -X local -i 30 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 1 ]
    [ "${lines[7]}" = "Exception: Unable to find the following arguments in the command: 'pos 30'" ]
}

@test "Testing good command" {
    run ./batch_launcher.py -C 10 -X local -i 2 -- grep "^>" test_data/HOT_100_reads.fasta > test_data/HOT_100_reads.local
    [ "$status" = 0 ]
    run diff <(sort test_data/HOT_100_reads.local) <(grep "^>" test_data/HOT_100_reads.fasta | sort)
    [ "$output" = "" ]
}




