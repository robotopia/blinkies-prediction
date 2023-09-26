# Predicting the orbital phase of GPM J1815-14

The GPM survey, combined with a Phase Dispersion Minimisation analysis, generates the following prediction for the orbital period:

![GPMJ1815-14_dutycycle.jpg](GPMJ1815-14_dutycycle.jpg)

Ordinarily, one would use this prediction to extrapolate to when the source should be on now.
However, in this case, the most constaining datum for when to observe is from the recent non-detection of the source from MeerKAT, which was observed between `2023-08-20T14:27:36.1` and `2023-08-20T19:39:02.4` UTC:
```
DATE-OBS= '2023-08-20T14:27:36.1'
DATE-OBS= '2023-08-20T14:49:36.4'
DATE-OBS= '2023-08-20T17:26:05.3'
DATE-OBS= '2023-08-20T17:48:11.3'
DATE-OBS= '2023-08-20T18:10:18.3'
DATE-OBS= '2023-08-20T18:32:27.8'
DATE-OBS= '2023-08-20T18:54:40.3'
DATE-OBS= '2023-08-20T19:16:52.9'
DATE-OBS= '2023-08-20T19:39:02.4'
```

If we assume this roughly 5-hour stretch falls neatly into the part of the alleged orbit when the source was off, and since the nominal off-period lasts `9.8*(1 - 0.53) = 4.6` hours, we can identify the middle of the observations with the middle of an off-period.
Predicting future ingress and egress times is then just a counting game: adding multiples of the estimated orbital period to the original assumed ingress and egress epochs.

The above logic is implemented in `J1815-14_predict.py`, which can be run with:
```
python J1815-14_predict.py
```

This script does not use command line options.
For now, the input parameters can be changed/tweaked directly in the script.
