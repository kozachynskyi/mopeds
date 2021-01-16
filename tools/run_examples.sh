#!/bin/bash

for d in $(find ./examples -name "*.py") ; do
    ignore=0
    for file in "./examples/dae_ode/hyfo_dae.py"
        do
            if [ "$file" == "$d" ]; then
                ignore=1
                break
            fi
        done

    echo "$d"
    if [[ $ignore = 0 ]]; then
        python "$d" >/dev/null 2>&1

        if [[ $? = 0 ]]; then
            echo "success"
        else
            echo "failure: $?"
        fi
    else
        echo "ignored"
    fi
done
