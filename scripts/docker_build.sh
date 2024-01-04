#!/bin/bash

# Check if Docker is running
if ! docker info >/dev/null 2>&1; then
    echo "Docker does not seem to be running, run it first and retry"
    exit 1
fi

# Build Docker images for all scripts
for dockerfile in $(find . -name 'Dockerfile_*'); do
    image_name=$(basename $dockerfile | cut -d'_' -f 2-)
    echo "Building image for $image_name"
    docker build -t $image_name -f $dockerfile .
    if [ $? -ne 0 ]; then
        echo "Failed to build image for $image_name"
        exit 1
    fi
done

echo "All images built successfully"
