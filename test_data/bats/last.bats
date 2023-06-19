#!/usr/bin/env bats



setup() {
    if [ ! -e test_data/.tst.genes.faa.prj ]; then 
        lastdb -p test_data/.tst.genes.faa test_data/genes.faa
    fi
}

@test "Test lastal" {
    run ./batch_launcher.py -N 10 -X local -- lastal -f blasttab test_data/.tst.genes.faa test_data/genes.faa
    [ "$status" = 0 ]
    run diff <(grep -v "^#" test_data/.tst.genes.faa.single | sort) <(grep -v "^#" test_data/genes.v.self.batch | sort)
    [ "$status" = 0 ]
}
