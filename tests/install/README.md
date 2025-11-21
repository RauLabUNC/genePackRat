### Intro
This folder contains the scripts necessary to produce an environment to test the install of `locusPackRat`
### How to use
First, build a `conda` environment from the `.yml` file. Then activate it:
```bash
conda env create -f tests/install/blank_slate_build.yml 
conda activate packrat_build_test
```