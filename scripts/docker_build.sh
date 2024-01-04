#!/bin/bash

# Check if docker is installed
if ! command -v docker &> /dev/null
then
    echo "Docker could not be found"
    exit
fi

# Build the docker image
docker build -t my-docker-image .

# Check if the build was successful
if [ $? -ne 0 ]
then
    echo "Docker build failed"
    exit 1
fi

echo "Docker build successful"
exit 0
```

File: tests/test_docker_build.sh

```bash
#!/bin/bash

# Run the docker_build.sh script
bash ../scripts/docker_build.sh

# Check if the script completed without errors
if [ $? -ne 0 ]
then
    echo "Test failed: docker_build.sh exited with errors"
    exit 1
fi

echo "Test passed: docker_build.sh completed without errors"
exit 0
