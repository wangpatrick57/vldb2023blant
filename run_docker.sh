#!/bin/bash
docker build -t blant .
docker run -v `pwd`:/home -w /home -it blant bash
