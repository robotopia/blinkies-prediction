# Predicting the orbital phase of GPM J1646-44

The orbital ephemeris (from Andy Wang) for this source is
```
T0 (MJD) 59391.3567 ± 0.0004
Pb (h) 5.267104 ± 3 × 10^−6
Ingress phase based on ASKAP 0.0502 ± 0.0016
Egress phase based on ASKAP 0.4497 ± 0.0016
```

The prediction of future ingress and egress times is then straightforward: count multiples of the orbital period from the reference ingress and egresses.
This is implemented in `J1646-44_predict.py`.
To run, type
```
python J1646-44_predict.py
```

This script does not use command line options.
For now, the input parameters can be changed/tweaked directly in the script.
