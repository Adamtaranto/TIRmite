# Code updates 

1. Re-structure package and split functions into modules
2. Use bedtools index to store seq on disk
3. Add docstrings 
4. Add logging
5. Add progress reporting
6. Update check for existing files methods
7. Change entry points to use sub-commands

# Dev tests
- Set up devcontainer for codespaces
- Create basic function tests and data
- Set up pytest action on push
- Test with Py 3.10

# Packaging 
- Add Dockerfile for package release (maybe diff to dev)
- Set up ghcr packaging actions on push tag
- Set up action to create pypi package on push tag

# Documentation
- Create basic jekyll page on gh-pages branch
