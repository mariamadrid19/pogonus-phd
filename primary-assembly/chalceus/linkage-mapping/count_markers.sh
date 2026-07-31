#!/bin/bash

printf "LG\tLOD10\tLOD11\tLOD12\tLOD15\n"

for lg in $(seq 1 12); do
    printf "%s" "$lg"

    for lod in 10 11 12 15; do
        count=$(awk -v lg="$lg" \
            '!/^#/ && $1 == lg {count++} END {print count+0}' \
            "map${lod}.txt")

        printf "\t%s" "$count"
    done

    printf "\n"
done
