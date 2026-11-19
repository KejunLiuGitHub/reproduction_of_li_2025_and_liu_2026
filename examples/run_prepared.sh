#!/usr/bin/env bash
# Run one newly prepared directory. Configure CP2K/MPI for your installation.
set -euo pipefail
: "${CP2K_BINARY:?Set CP2K_BINARY to the CP2K 9.1 executable}"
MPI_LAUNCHER=${MPI_LAUNCHER:-mpirun}
MPI_NP_FLAG=${MPI_NP_FLAG:--np}
MPI_TASKS=${MPI_TASKS:-64}
PYTHON=${PYTHON:-python3}
export OMP_NUM_THREADS=1
if [[ "$MPI_TASKS" != 64 ]]; then
  echo "These inputs were prepared for 64 MPI tasks." >&2
  exit 2
fi
if [[ $# != 1 ]]; then
  echo "Usage: run_prepared.sh NEW_PREPARED_DIRECTORY" >&2
  exit 2
fi
cd -- "$1"
if [[ -e STARTED_AT.txt || -e RETURN_CODE.txt ]]; then
  echo "This directory has already been launched; prepare a new one." >&2
  exit 2
fi
input=$("$PYTHON" -c 'import json; print(json.load(open("PREPARATION.json"))["input"])')
project=$("$PYTHON" -c 'import json; print(json.load(open("PREPARATION.json"))["project"])')
"$CP2K_BINARY" --version > ENGINE_VERSION.txt 2>&1
date -u '+%Y-%m-%dT%H:%M:%SZ' > STARTED_AT.txt
set +e
"$MPI_LAUNCHER" "$MPI_NP_FLAG" "$MPI_TASKS" "$CP2K_BINARY" -i "$input" -o "$project.out" > runner.stdout.txt 2> runner.stderr.txt
code=$?
set -e
printf '%s\n' "$code" > RETURN_CODE.txt
date -u '+%Y-%m-%dT%H:%M:%SZ' > ENDED_AT.txt
exit "$code"
