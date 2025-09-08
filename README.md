#  PHIDRA ![Version](https://img.shields.io/badge/version-2.0.0_stable-green)
**P**rotein **H**omology **I**dentification via **D**omain-**R**elated **A**rchitecture

A simple way to search and validate identifed Pfam domains of interest against a curated InterProScan Domain Architecture (IDA) file to check whether or not your proteins match a domain composition found in the Interpro Database.

## Description

A Python-based package of scripts that performs an initial homology search using `MMseqs2` to identify targeted proteins of interest in small or large datasets. Top hits are searched against the Pfam database using `pfam_scan`, and verified domains are checked and compared against a custom input InterProScan Domain Architecture (IDA) file relative to your target protein of interest. 
A recursive search is then performed using the full length proteins with validated IDAs as the subject database and the original input as the query, with initial matches filtered out. This process captures potentially more distant proteins that may have been missed in the initial homology search but are functionally relevant.

## Status

This `main` branch contains the redesigned codebase with **breaking changes**. 
Full documentation and usage examples are coming soon.

## Looking for Version 1?

- Browse the stable **[v1.x branch](https://github.com/zschreib/phidra/tree/v1.x)**  
- Or check the **[v1.0.0 release](https://github.com/zschreib/phidra/releases/tag/v1.0.0)**  
