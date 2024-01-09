#!/bin/bash

# Function to remove trailing period from the final atomtype column in a pdbqt file
remove_trailing_period() {
    input_pdbqt_file="$1"
    output_pdbqt_file="${input_pdbqt_file%.pdbqt}_fixed.pdbqt"

    sed 's/\(\S\+\)\.\s*$/\1/g' "$input_pdbqt_file" > "$output_pdbqt_file"
}

# Main function
main() {
    input_pdbqt_file="$1"

    # Call the remove_trailing_period function
    remove_trailing_period "$input_pdbqt_file"
}

# Execute the main function if the script is run directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
