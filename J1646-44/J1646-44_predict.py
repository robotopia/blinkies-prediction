import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo
from astroplan import Observer
from astropy.coordinates import EarthLocation
from astroplan import AltitudeConstraint
from astroplan import is_observable, is_always_observable, months_observable
from astroplan import FixedTarget
from astroplan.plots import plot_altitude
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from pytz import timezone
import matplotlib.pyplot as plt


gmrt_loc = EarthLocation.from_geodetic(lat=19.096517*u.deg, lon=74.049742*u.deg, height=650*u.m)
gmrt = Observer(name='uGMRT', location=gmrt_loc, timezone='Asia/Kolkata')
coords = SkyCoord("16 46 22.7", "-44 05 41", frame="fk5", unit=(u.hour, u.deg))
J1646 = FixedTarget(coords)

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#

# Altitude constraints at the uGMRT
ac = AltitudeConstraint(min=15*u.deg
)
# Change the range of desired MJD prediction here
MJD_start = Time(60250, format='mjd')
MJD_stop = Time(60260, format='mjd')

# Set the reference epoch (MJD)
T0 = Time(59391.3567, format='mjd')

# Output timezone
output_timezone = TimezoneInfo(utc_offset=5.5*u.hour)

# Offset reported times by some custom amount
# E.g. to start a one-hour dwell on-target, with the ingress or egress centred in the dwell, set to -0.5 hours
# Set to None or 0.0 for no offset
offset = -0.5 * u.hour

# Pre-determined orbital properties
Pb = 5.267104 * u.hour
ingress_epoch = T0 + 0.0502*Pb
egress_epoch = T0 + 0.4497*Pb

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# Calculate the number of orbital periods since the assumed egress/ingress epochs that fall within the given MJD range
egress_rotation_start = int(np.ceil((MJD_start - egress_epoch)/Pb))
egress_rotation_stop = int(np.ceil((MJD_stop - egress_epoch)/Pb))
ingress_rotation_start = int(np.ceil((MJD_start - ingress_epoch)/Pb))
ingress_rotation_stop = int(np.ceil((MJD_stop - ingress_epoch)/Pb))

# Calculate the actual ingress/egress times
egresses = np.arange(egress_rotation_start, egress_rotation_stop)*Pb + egress_epoch
ingresses = np.arange(ingress_rotation_start, ingress_rotation_stop)*Pb + ingress_epoch

# Apply the offset
if offset is not None:
    egresses += offset
    ingresses += offset

# Convert to specified timezone and print out
print(f"Egresses {offset.to('hour'):+f}:")
for egress in egresses:
    if is_observable(ac, gmrt, J1646, egress):
# Also needs to be observable one hour later
        endtime = egress + 1*u.hour
        if is_observable(ac, gmrt, J1646, endtime):
            print(egress.to_datetime(timezone=output_timezone))

print(f"Ingresses {offset.to('hour'):+f}:")
for ingress in ingresses:
    if is_observable(ac, gmrt, J1646, ingress):
# Also needs to be observable one hour later
        endtime = ingress + 1*u.hour
        if is_observable(ac, gmrt, J1646, endtime):
            print(ingress.to_datetime(timezone=output_timezone))
