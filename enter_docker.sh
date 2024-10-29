#!/bin/sh

docker run -it -v $(pwd):/workspace:rw mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04 /workspaces/
