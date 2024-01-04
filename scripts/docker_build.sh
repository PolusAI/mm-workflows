#!/bin/bash

# Define the Docker image name
DOCKER_IMAGE_NAME="mm-workflows"

# Check if the Docker image already exists
if [ "$(docker images -q ${DOCKER_IMAGE_NAME} 2> /dev/null)" == "" ]; then
    # If the Docker image does not exist, build it
    docker build -t ${DOCKER_IMAGE_NAME} .
else
    # If the Docker image already exists, print a message
    echo "Docker image '${DOCKER_IMAGE_NAME}' already exists."
fi
