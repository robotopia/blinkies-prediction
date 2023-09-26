import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#

# Change the range of desired MJD prediction here
MJD_start = Time(60250, format='mjd')
MJD_stop = Time(60260, format='mjd')

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
MeerKAT_start = Time('2023-08-20T14:27:36.1', format='isot', scale='utc')
MeerKAT_stop = Time('2023-08-20T19:39:02.4', format='isot', scale='utc') + obs_duration

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# How long is the source off?
off_time = Pb*(1 - duty_cycle)

# Calculate the middle of the non-detection time span
non_detection_ctr = MeerKAT_start + 0.5*(MeerKAT_stop - MeerKAT_start)
egress_epoch = non_detection_ctr - 0.5*off_time
ingress_epoch = non_detection_ctr + 0.5*off_time

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
print("Egresses:")
for egress in egresses:
    print(egress.to_datetime(timezone=output_timezone))

print("Ingresses:")
for ingress in ingresses:
    print(ingress.to_datetime(timezone=output_timezone))

