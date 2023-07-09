## Model curation

Scripts for curating different parts of the model.

Curations between each model release are gathered in consolidated curation scripts that are named after the model version that they generate. E.g. v8_6_1.m includes updates that regenerates model version 8.6.1. These scripts can replicate the curation process, but earlier model releases are best downloaded here: https://github.com/SysBioChalmers/yeast-GEM/releases.

TEMPLATEcuration.m is the template script that should be used to generate a new consolidated curation script after a new model has just been released.

In addition, this folder contains generic functions that are useful for model curation. These functions should be written in a generic format, so that they could be reused for future curations.
