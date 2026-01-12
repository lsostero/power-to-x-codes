Repository Description

This repository contains a set of MATLAB models for the techno-economic analysis (TEA) of Power-to-X pathways, with a focus on Power-to-Ammonia and Power-to-Methanol systems.
The code describes the full value chain, from renewable electricity generation to e-fuel production and final use in internal combustion engines.

Code Structure and Scope
1. E-fuel synthesis models

Core models for techno-economic analysis of synthetic fuel production:
ammoniasynthesis.m:Techno-economic model of green ammonia synthesis (Power-to-Ammonia).
methanolsythesis.m:Techno-economic model of methanol synthesis (Power-to-Methanol).

2. Fully renewable scenarios (wind-based systems)

Development of 100% renewable scenarios based on wind power as the primary energy source, including different storage and system integration strategies:
methanolsywindfarm.m: Wind-powered e-methanol production with system sizing and techno-economic evaluation.
greenammonia_windplushydrogenstorage.m: Wind-based green ammonia production including hydrogen storage to handle intermittency.

battery.m

Advanced green ammonia scenario with:
wind power
electrolysis
hydrogen storage
battery park
Includes integrated techno-economic analysis of the full system.

3. End-use: internal combustion engines fueled by e-fuels

Techno-economic and performance analysis of internal combustion engines (ICEs) operating on synthetic fuels:
methanolengine1.m: Techno-economic analysis of a methanol-fueled internal combustion engine.
ammoniaengine.m: Techno-economic analysis of an ammonia-fueled internal combustion engine.

Purpose
The repository enables:
end-to-end techno-economic assessment of Power-to-Ammonia and Power-to-Methanol pathways;
evaluation of wind-based, fully renewable e-fuel production systems;
analysis of the impact of hydrogen storage and batteries on system cost and performance;
linkage between fuel production TEA and final energy conversion in ICEs.
