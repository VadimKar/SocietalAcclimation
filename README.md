# SocietalAcclimation
Data and code to reproduce core manuscript analyses in Karatayev et al 2027 PRSB

Early attention speeds restoration: The dangers of acclimating to degraded ecosystems 

Authors: Vadim A. Karatayev, Robyn S. Wilson, D.G. Webster, Mark Axelrod, Chris T. Bauch, Madhur Anand

Published in Proceedings of the Royal Society B: Biological Sciences DOI 10.1098/rspb.


Distributed under an MIT license


Contents:

ErieEnvTrends.csv:  Data on trends in environmental variables in western Lake Erie.

IssueAttnCycles_Fig1.csv:  Data on trends in newspaper publications.

Windows 0.6 LITE.R:  Simulation and data analysis code



Software used:

R version 4.3.2 (2023-10-31)

R packages used:  pals v1.8,  rio v1.1.0,  fields v15.2,  viridisLite v0.4.2,  spam v2.10-0,  mgcv v1.9-0,  nlme v3.1-163,  nloptr v2.0.3,  abind v1.4-5,  deSolve v1.40     



ErieEnvTrends.csv: 

NonpointTP, UnmonitoredAdjTP, PointTP, AtmosphericTP, HuronTP, UndeterminedTP: annual total phosphorous loads into Lake Erie in tons separated by source. Data from ref 44.

AnnualDRPload_t: trends in soluble reactive phosphorous in western Lake Erie. Data from ref 44.

maxBloom_km2: Maximum algae bloom size (max across all blooms in a year). Data from Figure 1 of ref 39.

avgNdays, medNdays, avgSevere, medSevere: annual median and mean duration and severity levels of advisories posted for 68 Lake Erie beaches in Ohio in 2008-2020 from ref 42 (these statistics summarize across the 68 beaches).



IssueAttnCycles_Fig1.csv:

Area, Term, Year: search terms used for Newspapers.com search

URL: URL used to obtain results for our study, combining search terms

Hits: number of articles found on newspapers.com at that URL

URL_2026: corresponding URL on Newspapers.com as of November 12, 2026. Note that hit numbers no longer match our data, and in every case are higher, likely reflecting increased archive coverage since we ran our analysis.

YrSinceMax: Years since the maximum value of hits in the date range examined (maximum is at 0).

