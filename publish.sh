#!/bin/bash
cd "${1:-.}" || exit
echo $PWD
#maturin build --interpreter 3.9 --release --strip
#maturin publish --interpreter 3.9 --skip-existing
