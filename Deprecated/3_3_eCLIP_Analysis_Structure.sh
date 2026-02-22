#!/bin/bash

for fasta_file in *.fa; do
    # Extract CELL and PROTEIN names from the file name
    CELL=$(echo "$fasta_file" | cut -d'_' -f1)
    remaining=$(echo "$fasta_file" | cut -d'_' -f2-)
    if [[ "$remaining" =~ ([^_]+_[^_]+) ]]; then
        PROTEIN="${BASH_REMATCH[1]}"
    else
        PROTEIN=$(echo "$remaining" | cut -d'_' -f1)
    fi

    # Create folder if it doesn't exist
    folder_name="${CELL}_${PROTEIN}"
    mkdir -p "$folder_name"

    # Run RNAplfold
    # RNAplfold -W 25 -u 1 < "$fasta_file"
    # rm *.ps
    RNAfold -o < "$fasta_file"
    rm *.ps

    # Process output and move .bpp file to folder
    # for ps_file in *.ps; do
    #     sed -n '/%start of base pair probability data/,/showpage/p' "$ps_file" | awk '/%start of base pair probability data/{flag=1;next} /showpage/{flag=0} flag' > "$(basename "$ps_file" .ps).bpp"
    #     rm "$ps_file"
    #     mv "$(basename "$ps_file" .ps).bpp" "$folder_name/"
    # done

    # move to folder
    # for lunp_file in *_lunp; do
    #     mv "$lunp_file" "$folder_name/"
    # done

    for fold_file in *.fold; do
        mv "$fold_file" "$folder_name/"
    done

done

