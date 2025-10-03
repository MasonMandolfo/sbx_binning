-----------------------------------------------------------------

<img src="https://github.com/sunbeam-labs/sunbeam/blob/main/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_binning

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_binning/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_binning/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_binning)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_binning/)
<!-- Badges end -->

## Introduction

sbx_binning is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for binning assembled contigs from sbx_assembly into Metagenome-assembled genomes (MAGS) This pipeline uses binners metabat2 and VAMB and refines these bins using MAGScoT and CheckM2 for quality control. Before using the tool, you must install MAGScoT from github and configure the path/to/MAGScoT.

## Config

  - all_binning: runs the full binning pipeline
    
## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).
