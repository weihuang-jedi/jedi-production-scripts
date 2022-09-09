#!/bin/bash
#…or create a new repository on the command line
echo "# jedi-production-scripts" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M develop
git remote add origin https://github.com/weihuang-jedi/jedi-production-scripts.git
git push -u origin develop

#…or push an existing repository from the command line
#git remote add origin https://github.com/weihuang-jedi/jedi-production-scripts.git
#git branch -M develop
#git push -u origin develop

