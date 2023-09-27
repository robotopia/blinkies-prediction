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
coords = SkyCoord("18 15 56.8", "-14 16 34", frame="fk5", unit=(u.hour, u.deg))
J1815 = FixedTarget(coords)

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#

# Altitude constraints at the uGMRT
ac = AltitudeConstraint(min=30*u.deg)

# Change the range of desired MJD prediction here
MJD_start = Time(60240, format='mjd', location=gmrt_loc)
# We don't have much faith that the period is accurate so let's try to observe sooner rather than later
MJD_stop = Time(60240+100, format='mjd', location=gmrt_loc)

# Output timezone
output_timezone = TimezoneInfo(utc_offset=5.5*u.hour)

# Offset reported times by some custom amount (e.g. half an hour before the ingress/egress times)?
# Set to None or 0.0 for no offset
offset = -0.5 * u.hour

# Measurements from previous phase dispersion minimisation applied to MWA detections
Pb = 9.8 * u.hour
duty_cycle = 0.53

# The MeerKAT non-detection(s)
obs_duration = 550 * u.second
MeerKAT_start = Time('2023-08-20T14:27:36.1', format='isot', scale='utc', location=gmrt_loc)
MeerKAT_stop = Time('2023-08-20T19:39:02.4', format='isot', scale='utc', location=gmrt_loc) + obs_duration

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# How long is the source off?
off_time = Pb*(1 - duty_cycle)

# Calculate the middle of the non-detection time span
non_detection_ctr = MeerKAT_start + 0.5*(MeerKAT_stop - MeerKAT_start)
ingress_epoch = non_detection_ctr - 0.5*off_time # "ingress" means "going into eclipse"
egress_epoch = non_detection_ctr + 0.5*off_time # "egress" means "coming out of eclipse"

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

# For this source we only want to try to get some "on" measurements
# Convert to specified timezone and print out
#print(f"Egresses {offset.to('hour'):+f}:")
#for egress in egresses:
#    print(egress.to_datetime(timezone=output_timezone))

# Convert to specified timezone and print out
print(f"Egresses {offset.to('hour'):+f}:")
for egress in egresses:
    if is_observable(ac, gmrt, J1815, egress):
# Needs to start very close to (within five minutes of) a whole number of hours in IST
        t = egress.to_datetime(timezone=output_timezone)
        minutes = t.minute
        if minutes > 55 or minutes < 5:
    # Also needs to be observable five hours later
            endtime = egress + 5*u.hour
            if is_observable(ac, gmrt, J1815, endtime):
                 if minutes > 55:
                     hr = t.hour + 1
                 elif minutes < 5:
                     hr = t.hour
                 print("IST: {0}; LST: {1:2.3f}; Scheduling block: {2:4d}-{3:02d}-{4:02d} {5:02d}:00:00".format(t, egress.sidereal_time('apparent'), t.year, t.month, t.day, hr))

#print(f"Ingresses {offset.to('hour'):+f}:")
#for ingress in ingresses:
#    print(ingress.to_datetime(timezone=output_timezone))

