---
title: 'Solving the "many variables" problem in MICE with principal component regression'
authors:
  - admin
  - Kyle M. Lang
  - Klaas Sijtsma
  - Tim Reeskens
author_notes:
  - 'main author'
  - 'main advisor'
  - 'advisor'
  - 'advisor'
date: '2022-06-30T00:10:00Z'
doi: 'https://doi.org/10.48550/arXiv.2206.15107'
publishDate: '2017-01-01T00:00:00Z'
publication_types: ['3']
publication: arXiv
publication_short: 
abstract: Multiple Imputation (MI) is one of the most popular approaches to addressing missing values in questionnaires and surveys. MI with multivariate imputation by chained equations (MICE) allows flexible imputation of many types of data. In MICE, for each variable under imputation, the imputer needs to specify which variables should act as predictors in the imputation model. The selection of these predictors is a difficult, but fundamental, step in the MI procedure, especially when there are many variables in a data set. In this project, we explore the use of principal component regression (PCR) as a univariate imputation method in the MICE algorithm to automatically address the "many variables" problem that arises when imputing large social science data. We compare different implementations of PCR-based MICE with a correlation-thresholding strategy by means of a Monte Carlo simulation study and a case study. We find the use of PCR on a variable-by-variable basis to perform best and that it can perform closely to expertly designed imputation procedures.
summary: Comparison of different implementations of PCR-based MICE
tags: ["missing values", "imputation"]
featured: no
url_pdf: 'https://arxiv.org/pdf/2206.15107.pdf'
url_dataset: 'https://github.com/amices/mice/blob/master/data/fdd.rda'
url_poster: ''
url_project: ''
url_slides: ''
url_source: ''
url_video: ''
links:
- name: Code - simulation study
  url: https://github.com/EdoardoCostantini/mipca_compare
- links:
  name: Code - sase study
  url: https://github.com/EdoardoCostantini/fireworks
projects:
  - phd-research
---

### Results dashboard

Here you can explore the results reported in the simulation study. 
Thanks to this shiny app you can actually check all the additional results that couldn't be included in the main text.

<iframe height="900" width="100%" frameborder="no" src="https://edoardocostantini.shinyapps.io/shiny-mi-pca-plot/?_ga=2.187725984.929064409.1659041054-1213691852.1658930327"> </iframe>