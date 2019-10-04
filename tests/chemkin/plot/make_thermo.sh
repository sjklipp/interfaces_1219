#!/bin/bash

rm thm_plots/*
python test__thermo.py
evince thm_plots/all.pdf

