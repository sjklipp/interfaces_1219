#!/bin/bash

rm rate_plots/*
python test__rates.py
evince rates_plots/all.pdf

