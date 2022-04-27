---
draft: true
title: "{{ replace .Name "-" " " | title }}"
subtitle: ""
summary: ""
authors: ["admin"]
tags: []
categories: []
date: {{ .Date }}
lastmod: {{ .Date }}
featured: false
image:
  caption: ""
  focal_point: ""
  preview_only: false
projects: []
output:
  blogdown::html_page:
    toc: true
    toc_depth: 4
    number_sections: true
---
