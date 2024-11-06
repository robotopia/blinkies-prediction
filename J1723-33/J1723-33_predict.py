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
coords = SkyCoord("17 23 38.57", "-33 41 57.34", frame="fk5", unit=(u.hour, u.deg))
J1723 = FixedTarget(coords)

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#

# Altitude constraints at the uGMRT
ac = AltitudeConstraint(min=30*u.deg)

# Change the range of desired MJD prediction here
MJD_start = Time("2024-11-01T00:00:00", format='isot', scale="utc", location=gmrt_loc)
# We don't have much faith that the period is accurate so let's try to observe sooner rather than later
MJD_stop = Time("2025-01-31T00:00:00", format='isot', scale="utc", location=gmrt_loc)

# Output timezone
output_timezone = TimezoneInfo(utc_offset=5.5*u.hour)

# Offset reported times by some custom amount (e.g. half an hour before the ingress/egress times)?
# Set to None or 0.0 for no offset
offset = -0.5 * u.hour

#T0 = Time(59805.48950, format='mjd', scale='utc', location=gmrt_loc)# start of eclipse
T0 = Time(59789.56601, format='mjd', scale='utc', location=gmrt_loc)# middle of eclipse

# Measurements from previous phase dispersion minimisation applied to MWA detections
Pb = 0.9451 * u.day

# Time that it spends > 50% 'on'
duty_cycle = 0.65

# Duration of the ingresses/egresses -- pad it a little bit to make sure we catch a full one
duration = 0.15*Pb

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# How long between 50% on?
off_time = Pb*(1 - duty_cycle)

# Noting that middle of ingress is actually earlier than Flora's T0 by about 0.1 orbital phase
# So middle of ingress is current T0 plus half of the duty_cycle; beginning is minus half of the duration
ingress_epoch = T0 + 0.5*duty_cycle*Pb - 0.5*duration
# And middle of egress is current T0 minus half of the duty_cyle; beginning is minus half of the duration
egress_epoch = T0 - 0.5*duty_cycle*Pb - 0.5*duration

# Calculate the number of orbital periods since the assumed egress/ingress epochs that fall within the given MJD range
egress_rotation_start = int(np.ceil((MJD_start - egress_epoch)/Pb))
egress_rotation_stop = int(np.ceil((MJD_stop - egress_epoch)/Pb))
ingress_rotation_start = int(np.ceil((MJD_start - ingress_epoch)/Pb))
ingress_rotation_stop = int(np.ceil((MJD_stop - ingress_epoch)/Pb))

# Calculate the actual ingress/egress times
egresses = np.arange(egress_rotation_start, egress_rotation_stop)*Pb + egress_epoch
ingresses = np.arange(ingress_rotation_start, ingress_rotation_stop)*Pb + ingress_epoch

# Apply the offset (to allow time for phase calibration etc)
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
    if is_observable(ac, gmrt, J1723, egress):
# Will have to relax the timing constraint and then the offset because we're not finding anything
        t = egress.to_datetime(timezone=output_timezone)
        minutes = t.minute
        if minutes > 30 or minutes < 30:
    # Also needs to be observable "duration" hours later
            starttime = egress
            endtime = egress + duration.to(u.hour)
            t = starttime.to_datetime(timezone=output_timezone)
            if is_observable(ac, gmrt, J1723, starttime) and is_observable(ac, gmrt, J1723, endtime):
                 if minutes > 30:
                     hr = t.hour + 1
                 elif minutes < 30:
                     hr = t.hour
                 print("IST: {0}; LST: {1:2.3f}; Scheduling block: {2:4d}-{3:02d}-{4:02d} {5:02d}:00:00".format(t, starttime.sidereal_time('apparent'), t.year, t.month, t.day, hr))

print(f"Ingresses {offset.to('hour'):+f}:")
for ingress in ingresses:
    if is_observable(ac, gmrt, J1723, ingress):
# Needs to start very close to (within ten minutes of) a whole number of hours in IST
        t = ingress.to_datetime(timezone=output_timezone)
        minutes = t.minute
        if minutes > 30 or minutes < 30:
    # Also needs to be observable "duration" hours later
            starttime = ingress
            endtime = ingress + duration.to(u.hour)
            t = starttime.to_datetime(timezone=output_timezone)
            if is_observable(ac, gmrt, J1723, starttime) and is_observable(ac, gmrt, J1723, endtime):
                 if minutes > 30:
                     hr = t.hour + 1
                 elif minutes < 30:
                     hr = t.hour
                 print("IST: {0}; LST: {1:2.3f}; Scheduling block: {2:4d}-{3:02d}-{4:02d} {5:02d}:00:00".format(t, starttime.sidereal_time('apparent'), t.year, t.month, t.day, hr))

