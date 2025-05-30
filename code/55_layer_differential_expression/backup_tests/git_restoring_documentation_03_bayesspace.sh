#!/bin/bash

# Steps to revert a specific file to a specifi commit hash

## Get commit history and check possible disruptive change

git log -- 03_model_BayesSpace.R > gitlog_03_model_BayesSpace.txt
cat gitlog_03_model_BayesSpace.txt
git show 4165b394b464184ac9cf3bb309bf846777c7eccd


## Now, revert all commits to the specific desired point
## Note:
## Git reverts in reverse order internally, which avoids dependency issues between commits
## Git reverts in reverse order internally, which avoids dependency issues between commits

## To revert the file 03_model_BayesSpace.R to its state from commit f746b10

# That range includes 3 commits in order:
    
f746b10193e83a4df7c23b8ee2e2005f9bcfd279 # ← target state (Mon Jul 1 13:21:05 2024 -0400)

bf06284e351e9c062f24635dc978bbfd033ed8db 

57862992717038a376a7ab312c48e4d2907f52fd

12d19bb268a5bbdbd49e5d3d05f18151734718c6

65befd616e72cd212e9441b093a8faa7fdaf4d45 # ← most recent (Fri May 23 08:59:02 2025 -0400)

# To revert those 3 newer commits (assuming they’re linear)
# Each git revert creates a new commit that undoes just the effect of that single commit
# Always revert from newest to oldest, or else later reverts might fail due to conflicts from earlier changes already being undone.

git revert 65befd616e72cd212e9441b093a8faa7fdaf4d45

git revert 12d19bb268a5bbdbd49e5d3d05f18151734718c6
# Changes to be committed:
#       modified:   code/04_harmony_BayesSpace/05-multi_gene_selection.R -> make copy 
#       modified:   code/05_layer_differential_expression/03_model_BayesSpace.R -> ok
git restore --staged code/04_harmony_BayesSpace/05-multi_gene_selection.R
git checkout -- code/04_harmony_BayesSpace/05-multi_gene_selection.R

git revert 57862992717038a376a7ab312c48e4d2907f52fd

git revert bf06284e351e9c062f24635dc978bbfd033ed8db

git revert f746b10193e83a4df7c23b8ee2e2005f9bcfd279


## Alternative: if the commits are consecutive and in order
## If you're not sure that commits affects only one file, avoid this rever option and check one by one
    
git revert <oldest-commit>^..<newest-commit>
   
git add 03_model_BayesSpace.R

git commit -m "Revert 03_model_BayesSpace.R to version from commit f746b10 (July 1, 2024)"


## The script may now match its state at an older commit 
## To check if your script matches the exact file version from a past commit:

git diff f746b10 -- code/05_layer_differential_expression/03_model_BayesSpace.R

## no output, it means the current version of the file matches exactly what it was in f746b10

## Alternativately

git show bf06284:code/05_layer_differential_expression/03_model_BayesSpace.R > old_version.R
diff code/05_layer_differential_expression/03_model_BayesSpace.R old_version.R



