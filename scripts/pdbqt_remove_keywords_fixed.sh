#!/bin/bash

# Function to remove the MODEL and ENDMDL keywords from a PDBQT file
remove_keywords() {
    input_pdbqt_file="$1"
    output_pdbqt_file="${input_pdbqt_file%.pdbqt}_fixed.pdbqt"

    sed '/^MODEL\|^ENDMDL/d' "$input_pdbqt_file" > "$output_pdbqt_file"
}

# Main function
main() {
    input_pdbqt_file="$1"

    # Call the remove_keywords function
    remove_keywords "$input_pdbqt_file"
}

# Execute the main function if the script is run directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
