#!/bin/bash

# Function to remove the ROOT, ENDROOT, and TORSDOF keywords from a PDBQT file
remove_flex_keywords() {
    input_pdbqt_file="$1"
    output_pdbqt_file="${input_pdbqt_file%.pdbqt}_fixed.pdbqt"

    sed '/^ROOT\|^ENDROOT\|^TORSDOF/d' "$input_pdbqt_file" > "$output_pdbqt_file"
}

# Main function
main() {
    input_pdbqt_file="$1"

    # Call the remove_flex_keywords function
    remove_flex_keywords "$input_pdbqt_file"
}

# Execute the main function if the script is run directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
