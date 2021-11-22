# Planetary-Defence
Git repository for [Kelvins Planetary Defence Challenge](http://kelvins.esa.int/planetary-defence/)
=========================


This can also act as a light curve analysis code, with both feature extraction and regeression algorithms included in the package.
--------------------------


Code Description: 

1. _allfeatures.m_, _newfeaturelist1903.m_ and _newfrequencies2607.m_ are the primary feature extraction codes. They should be executed in reverse order and will give the final feature list to be used for the regression algorithm to predict the shape and orbital paramters of the binary asteroid before and after collision.

2. _kelvins.ipynb_ is a jupyter notebook in python which uses sklearn packages to implement a feature regression and prediction algorithm

3. _frequencyplomb.py_ is used to get the PSD of the lightcurve data from python's LombScargle implementation.

4. _reconstructed.m_ uses PSD of the light curve to reconstruct a cleaner simpler version of the light curves, with less high frequency components.

The lightcurves are included in the lightcurve package and can also be accessed from the competition portal.