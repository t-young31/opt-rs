#!/bin/bash -l

for minor_version in {7..10}; do
    conda create -n test "python=3.$minor_version" --yes
    conda activate test
    pip install maturin
    maturin build --release --strip --universal2
    conda deactivate
done
