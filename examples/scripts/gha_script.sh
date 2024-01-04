#!/bin/bash

function execute_command() {
    output=$(echo "CI status received from 'PolusAI/workflow-inference-compiler'")
    status=$?

    if [ $status -ne 0 ]; then
        echo "Error: Command execution failed with error: $output"
        exit $status
    else
        echo $output
        exit 0
    fi
}

execute_command
