#!/usr/bin/env bats

teardown() {
    rm -f "bats-retry-${BATS_TEST_NUMBER}.stderr" \
        "bats-retry-${BATS_TEST_NUMBER}.failures" \
        "bats-failure-${BATS_TEST_NUMBER}.stderr" \
        "bats-failure-${BATS_TEST_NUMBER}.failures"
}

@test "Testing missing chunk flag" {
    run ./batch_launcher.py -X local -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 2 ]
}

@test "Rejecting unsupported SGE option" {
    run ./batch_launcher.py -X sge -C 10 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 2 ]
    [[ "$output" == *"sge option is no longer supported"* ]]
}

@test "Rejecting unsupported SLURM option" {
    run ./batch_launcher.py -X slurm -C 10 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 2 ]
    [[ "$output" == *"slurm option is no longer supported"* ]]
}

@test "Testing missing input designation" {
    run ./batch_launcher.py -C 10 -X local -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 1 ]
}

@test "Testing wrong input designation" {
    run ./batch_launcher.py -C 10 -X local -i 30 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 1 ]
    [[ "$output" == *"Exception: Unable to find the following arguments in the command: 'pos 30'"* ]]
}

@test "Testing good command" {
    run ./batch_launcher.py -C 10 -i 2 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 0 ]
    [ "${#lines[@]}" = 100 ]
}

@test "Testing no splitting" {
    run ./batch_launcher.py -N 1 -X local -i 2 -- grep "^>" test_data/HOT_100_reads.fasta
    [ "$status" = 0 ]
    [ "${#lines[@]}" = 100 ]
}

@test "Reading chunked FASTA input from stdin" {
    run bash -c 'printf ">read-1\nAC\n>read-2\nGT\n" | ./batch_launcher.py -C 1 -T fasta -- grep "^>"'
    [ "$status" = 0 ]
    [ "${#lines[@]}" = 2 ]
    [ "${lines[0]}" = ">read-1" ]
    [ "${lines[1]}" = ">read-2" ]
}

@test "Mapping and merging an explicit output file" {
    input="$BATS_TEST_TMPDIR/input.fasta"
    output_path="$BATS_TEST_TMPDIR/merged.fasta"
    command="$BATS_TEST_TMPDIR/copy-flags"
    cat > "$input" <<'FASTA'
>read-1
AC
>read-2
GT
>read-3
AA
>read-4
TT
FASTA
    cat > "$command" <<'SH'
#!/bin/sh
while [ "$#" -gt 0 ]; do
    case "$1" in
        --input) input="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        *) shift ;;
    esac
done
cp "$input" "$output"
SH
    chmod +x "$command"

    run ./batch_launcher.py -N 2 -i=--input -o=--output -- "$command" --input "$input" --output "$output_path"
    [ "$status" = 0 ]
    diff "$input" "$output_path"
}

@test "Retrying a fragment that fails once" {
    input="$BATS_TEST_TMPDIR/input.fasta"
    marker="$BATS_TEST_TMPDIR/failed-once"
    command="$BATS_TEST_TMPDIR/retry-once"
    result="$BATS_TEST_TMPDIR/result.fasta"
    stderr_path="$BATS_TEST_TMPDIR/stderr"
    cat > "$input" <<'FASTA'
>read-1
AC
>read-2
GT
>read-3
AA
>read-4
TT
FASTA
    cat > "$command" <<'SH'
#!/bin/sh
if (set -C; : > "$BATCH_TEST_MARKER") 2>/dev/null; then
    exit 23
fi
cat "$1"
SH
    chmod +x "$command"

    run bash -c 'BATCH_TEST_MARKER="$1" ./batch_launcher.py -j "bats-retry-$6" -N 4 -r 2 -i 1 -- "$2" "$3" > "$4" 2> "$5"' \
        _ "$marker" "$command" "$input" "$result" "$stderr_path" "$BATS_TEST_NUMBER"
    [ "$status" = 0 ]
    [ -f "$marker" ]
    grep -q "Error code 23" "$stderr_path"
    sort "$input" > "$BATS_TEST_TMPDIR/expected.sorted"
    sort "$result" > "$BATS_TEST_TMPDIR/result.sorted"
    diff "$BATS_TEST_TMPDIR/expected.sorted" "$BATS_TEST_TMPDIR/result.sorted"
}

@test "Propagating failure when retries are disabled" {
    input="$BATS_TEST_TMPDIR/input.fasta"
    command="$BATS_TEST_TMPDIR/always-fail"
    cat > "$input" <<'FASTA'
>read-1
AC
>read-2
GT
FASTA
    cat > "$command" <<'SH'
#!/bin/sh
exit 23
SH
    chmod +x "$command"

    run ./batch_launcher.py -j "bats-failure-$BATS_TEST_NUMBER" -N 2 -r 0 -i 1 -- "$command" "$input"
    [ "$status" -ne 0 ]
    [[ "$output" == *"Error code 23"* ]]
}

@test "Respecting the local task throttle" {
    input="$BATS_TEST_TMPDIR/input.fasta"
    command="$BATS_TEST_TMPDIR/track-concurrency"
    state="$BATS_TEST_TMPDIR/concurrency"
    printf '0 0\n' > "$state"
    for i in {1..8}; do
        printf '>read-%s\nACGT\n' "$i"
    done > "$input"
    cat > "$command" <<'PY'
#!/usr/bin/env python3
import fcntl
import sys
import time

state_path = sys.argv[2]
with open(state_path, "r+") as state:
    fcntl.flock(state, fcntl.LOCK_EX)
    active, maximum = map(int, state.read().split())
    active += 1
    maximum = max(maximum, active)
    state.seek(0)
    state.truncate()
    state.write("{} {}\n".format(active, maximum))
    state.flush()

time.sleep(0.25)

with open(state_path, "r+") as state:
    fcntl.flock(state, fcntl.LOCK_EX)
    active, maximum = map(int, state.read().split())
    state.seek(0)
    state.truncate()
    state.write("{} {}\n".format(active - 1, maximum))
    state.flush()

with open(sys.argv[1]) as fragment:
    sys.stdout.write(fragment.read())
PY
    chmod +x "$command"

    run ./batch_launcher.py -N 8 -R 2 -i 1 -- "$command" "$input" "$state"
    [ "$status" = 0 ]
    [ "${#lines[@]}" = 16 ]
    read -r active maximum < "$state"
    [ "$active" = 0 ]
    [ "$maximum" -gt 1 ]
    [ "$maximum" -le 2 ]
}
