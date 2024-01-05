#!/bin/bash

# This is a hypothetical representation of the docker_build.sh script. The actual script may vary.

echo "Building Docker image..."

if [ "$1" == "" ]; then
    echo "No argument supplied"
    exit 1
fi

IMAGE_NAME=$1

docker build -t $IMAGE_NAME .

if [ $? -eq 0 ]; then
    echo "Docker image built successfully."
else
    echo "Failed to build Docker image."
    exit 1
fi
