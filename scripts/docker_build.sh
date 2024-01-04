#!/bin/bash

# Define the Docker image name
DOCKER_IMAGE_NAME="mm-workflows"

# Check if Docker is running
if ! docker info >/dev/null 2>&1; then
    echo "Docker does not seem to be running, run it first and retry"
    exit 1
fi

# Build the Docker image
docker build -t ${DOCKER_IMAGE_NAME} .

# Check if the Docker image was successfully built
if [[ "$(docker images -q ${DOCKER_IMAGE_NAME} 2> /dev/null)" == "" ]]; then
    echo "Error building Docker image. Please check the Dockerfile and try again."
    exit 1
fi

# If everything is successful, print a success message
echo "Docker image '${DOCKER_IMAGE_NAME}' built successfully."
