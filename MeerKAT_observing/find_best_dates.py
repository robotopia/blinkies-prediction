import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

from astroplan import Observer
from astroplan import AltitudeConstraint, PhaseConstraint
from astroplan import is_observable, is_always_observable, months_observable
from astroplan import FixedTarget
from astroplan.plots import plot_altitude
from astroplan import PeriodicEvent

from pytz import timezone
import matplotlib.pyplot as plt

mkt = Observer.at_site("salt")

coords = {"J1646": FixedTarget(SkyCoord("16 46 22.7", "-44 05 41", frame="fk5", unit=(u.hour, u.deg)), name="J1646"),
          "J1723" : FixedTarget(SkyCoord("17 23 38.5", "-33 42 00", frame="fk5", unit=(u.hour, u.deg)), name="J1723"),
          "J1728" : FixedTarget(SkyCoord("17 28 12.12", "-46 08 01.49", frame="fk5", unit=(u.hour, u.deg)), name="J1728"),
          "J1734" : FixedTarget(SkyCoord("17 34 33.9", "-28 07 07", frame="fk5", unit=(u.hour, u.deg)), name="J1734"),
          "J1740" : FixedTarget(SkyCoord("17 40 16.16", "-26 50 28.77", frame="fk5", unit=(u.hour, u.deg)), name="J1740"),
          "J1752" : FixedTarget(SkyCoord("17 52 26.2", "-30 33 43", frame="fk5", unit=(u.hour, u.deg)), name="J1752"),
          "J1815" : FixedTarget(SkyCoord("18 15 56.8", "-14 16 34", frame="fk5", unit=(u.hour, u.deg)), name="J1815")}

# T0 = middle of eclipse, in MJD
# THESE ARE PLACEHOLDERS UNTIL CONFIRMED BY FLORA
T0s = {"J1646" : Time(59700, format='mjd', scale='utc', location=mkt.location),
       "J1723" : Time(59788.54, format='mjd', scale='utc', location=mkt.location),
       "J1728" : Time(59700, format='mjd', scale='utc', location=mkt.location),
       "J1734" : Time(59765.611, format='mjd', scale='utc', location=mkt.location),
       "J1740" : Time(59700, format='mjd', scale='utc', location=mkt.location),
       "J1752" : Time(59791.54, format='mjd', scale='utc', location=mkt.location),
       "J1815" : Time(59732.73, format='mjd', scale='utc', location=mkt.location)}

# Orbital periods
# GOOD FOR FLORA TO CHECK THESE
Pbs = {"J1646" : 5.26703*u.hour,
       "J1723" : 22.6824*u.hour,
       "J1728" : 5.04983*u.hour,
       "J1734" : 10.93128*u.hour,
       "J1740" : 16.50614*u.hour,
       "J1752" : 17.8632*u.hour,
       "J1815" : 9.82272*u.hour}

# duty cycle of ON phase
duty_cycles = {"J1646" : 0.5,
       "J1723" : 0.8,
       "J1728" : 0.6,
       "J1734" : 0.7,
       "J1740" : 0.2,
       "J1752" : 0.4,
       "J1815" : 0.5}

# Astroplan Periodic Event objects
pe = {}
for obj in T0s:
    pe[obj] = PeriodicEvent(epoch=T0s[obj], period=Pbs[obj])

# Altitude constraints at MeerKAT
ac = AltitudeConstraint(min=15*u.deg)

# Observing constraints (right part of orbital phase, and above the horizon)
constraints = {}
for obj in T0s:
    constraints[obj] = [PhaseConstraint(pe[obj], min=0.5 - (duty_cycles[obj]/2), max=0.5 + (duty_cycles[obj]/2)), ac]

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#


# Change the range of desired MJD prediction here
MJD_start = Time("2024-11-15T00:00:00", format='isot', scale="utc", location=mkt.location)
#MJD_stop = Time("2025-11-15T00:00:00", format='isot', scale="utc", location=mkt.location)
# For debugging
MJD_stop = Time("2024-11-22T00:00:00", format='isot', scale="utc", location=mkt.location)

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# For each day, calculate whether all of the sources are observable or if there is a showstopper:
ok = []
for mjd in np.arange(MJD_start.mjd, MJD_stop.mjd):
    for obj in T0s:
        if not is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd, format='mjd', scale='utc'), Time(mjd+1, format='mjd', scale='utc')]):
            print(f"Can't observe {obj} on {mjd}")
    ok.append(mjd) 



# Remove any day where any source is completely in eclipse the whole observation





# Output timezone = UTC
#output_timezone = TimezoneInfo(utc_offset=0*u.hour)
#
## Offset reported times by some custom amount (e.g. half an hour before the ingress/egress times)?
## Set to None or 0.0 for no offset
#offset = 0.0 * u.hour
#
#
## How long between 50% on?
#off_time = Pb*(1 - duty_cycle)
#
## T0 is the middle of the eclipse
#
## So middle of ingress is current T0 minus half of 1-duty_cycle; beginning is minus half of the duration
#ingress_epoch = T0 - 0.5*(1-duty_cycle)*Pb - 0.5*duration
## And middle of egress is current T0 plus half of 1-duty_cyle; beginning is minus half of the duration
#egress_epoch = T0 + 0.5*(1-duty_cycle)*Pb - 0.5*duration
#
## Calculate the number of orbital periods since the assumed egress/ingress epochs that fall within the given MJD range
#egress_rotation_start = int(np.ceil((MJD_start - egress_epoch)/Pb))
#egress_rotation_stop = int(np.ceil((MJD_stop - egress_epoch)/Pb))
#ingress_rotation_start = int(np.ceil((MJD_start - ingress_epoch)/Pb))
#ingress_rotation_stop = int(np.ceil((MJD_stop - ingress_epoch)/Pb))
#
## Calculate the actual ingress/egress times
#egresses = np.arange(egress_rotation_start, egress_rotation_stop)*Pb + egress_epoch
#ingresses = np.arange(ingress_rotation_start, ingress_rotation_stop)*Pb + ingress_epoch
#
## Apply the offset (to allow time for phase calibration etc)
#if offset is not None:
#    egresses += offset
#    ingresses += offset
#
## For this source we only want to try to get some "on" measurements
## Convert to specified timezone and print out
##print(f"Egresses {offset.to('hour'):+f}:")
##for egress in egresses:
##    print(egress.to_datetime(timezone=output_timezone))
#
## How many minutes you are happy to start away from integer IST scheduling blocks
#flex = 20
#
## Convert to specified timezone and print out
#print(f"Egresses {offset.to('hour'):+f}:")
#for egress in egresses:
#    if is_observable(ac, gmrt, J1723, egress):
## Will have to relax the timing constraint and then the offset because we're not finding anything
#        t = egress.to_datetime(timezone=output_timezone)
#        minutes = t.minute
#        if minutes > 60-flex or minutes < flex:
#    # Also needs to be observable "duration" hours later
#            starttime = egress
#            endtime = egress + duration.to(u.hour)
#            t = starttime.to_datetime(timezone=output_timezone)
#            if is_observable(ac, gmrt, J1723, starttime) and is_observable(ac, gmrt, J1723, endtime):
#                 if minutes > 60-flex:
#                     hr = t.hour + 1
#                 elif minutes < flex:
#                     hr = t.hour
#                 print("IST: {0}; LST: {1:2.3f}; Scheduling block: {2:4d}-{3:02d}-{4:02d} {5:02d}:00:00".format(t, starttime.sidereal_time('apparent'), t.year, t.month, t.day, hr))
#
#print(f"Ingresses {offset.to('hour'):+f}:")
#for ingress in ingresses:
#    if is_observable(ac, gmrt, J1723, ingress):
## Needs to start very close to (within ten minutes of) a whole number of hours in IST
#        t = ingress.to_datetime(timezone=output_timezone)
#        minutes = t.minute
#        if minutes > 60-flex or minutes < flex:
#    # Also needs to be observable "duration" hours later
#            starttime = ingress
#            endtime = ingress + duration.to(u.hour)
#            t = starttime.to_datetime(timezone=output_timezone)
#            if is_observable(ac, gmrt, J1723, starttime) and is_observable(ac, gmrt, J1723, endtime):
#                 if minutes > 60-flex:
#                     hr = t.hour + 1
#                 elif minutes < flex:
#                     hr = t.hour
#                 print("IST: {0}; LST: {1:2.3f}; Scheduling block: {2:4d}-{3:02d}-{4:02d} {5:02d}:00:00".format(t, starttime.sidereal_time('apparent'), t.year, t.month, t.day, hr))

