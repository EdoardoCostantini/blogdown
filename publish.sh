#!/bin/bash
# Project:  blogdown
# Topic:    Bash script to build and publish the website to github
# Author:   Edoardo Costantini
# Created:  2022-04-07
# Modified: 2022-04-07

# Clean current website
rm -r ~/projects/blogdown/public/*
rm -r ~/projects/edoardocostantini.github.io/*

# Build website
R -e 'blogdown::build_site()'

# Copy public to github website location
cp -r ~/projects/blogdown/public/. ~/projects/edoardocostantini.github.io/

# Move to git folder for website location
cd ~/projects/edoardocostantini.github.io/

# Stage all changes
git add -A

# Commit changes
git commit -m "automatic website update"

# Push changes to GitHub
git push

# Return to previous location
cd ~/projects/blogdown