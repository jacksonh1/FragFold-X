#!/bin/bash
# 1. generate the input param yamls (no GPU needed)
python generate_params.py
# 2. run each generated prediction, resolving relative paths against this directory
for params in params/*.yaml; do
    fragfoldx --input_params "$params" --root .
done
