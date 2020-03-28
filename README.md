# Computational models to simulate epidemics

This repository contains MATLAB codes for ODE and individual-based stochastic versions of the classic SIR model. These models are for educational use and are not appropriate for healthcare decision-making.

## MATLAB code
SIRode.m:     ordinary differential equations defining the SIR model
SIRcompute.m: script that runs basic simulations with SIRode.m
SIRsensitivity.m: runs various sensitivity analyses- sensitivity coefficients, response curves, sensitivity matrix, and uncertainity analysis
SIRabm.m: stochastic individual-based SIR model (non-spatial)

saucerman_epidemicsLecture.pdf: a subset of slides I used for a lecture in BME 8315: Systems Biology and Multi-scale Models. The focus of the lecture was on parallels in approaches to modeling across scales, and on how the assumptions and goals underlying model development drive appropriate mathematical representation.

For a more visual introduction to epidemic modeling intended for elementary, middle, or high-school students, see our project Go Viral at:
https://sites.google.com/virginia.edu/goviral

This article has a very nice description of the development, fundamentals, and applications of the SIR model:
Weiss, Howard Howie. “The SIR Model and the Foundations of Public Health.” Materials Matematics, 2013, 0001–17. http://mat.uab.cat/matmat/PDFv2013/v2013n03.pdf
