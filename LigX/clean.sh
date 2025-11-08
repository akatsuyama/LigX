#!/bin/bash

find . -name ".DS_Store" | xargs rm
find . -name ".DS_Store" | xargs rm -r
dot_clean .
