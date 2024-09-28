# ecef_to_eci.py

# Usage: python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km
# Test Case: py ecef_to_eci.py 2054 4 29 11 29 03.3 6778.136999646678 -0.030015095972430572 3838.027968
# Converting from ecef to eci positions

# Parameters:
#  year : int 
#   year
#  month : int 
#   month
#  day : int 
#   day
#  hour : int 
#   hour
#  minute : int 
#   min
#  second : float
#   sec
#  ecef_x_km : int | float | str
#   ECEF X position in km
#  ecef_y_km : int | float | str
#   ECEF Y position in km
#  ecef_z_km : int | float | str
#   ECEF Z position in km


# Output:
#  eci : list
# eci position vector

# Written by Riley Parsons

import sys
import math

# "constants"
w = 7.292115e-5

# helper functions  
def ymdhms_tojd(y : int, m : int, d: int, hr: int, min: int, sec: float):
  jd = d - 32075 \
      + int(1461 * (y + 4800 + int((m - 14)/12))/4) \
      + int(367 * (m - 2 - int((m-14)/12) * 12)/12) \
      - int(3 * int(((y + 4900 + int((m-14)/12))/100))/4)
  
  jd_midnight = jd - 0.5
  d_frac = (sec + 60 * (min + 60 * hr ))/86400
  jd_frac = jd_midnight + d_frac
    
  return jd_frac

def jd_to_tut(jd_frac):
  return (jd_frac - 2451545.0)/36525

def tut_to_GMST_rads(tut):
  gmst_angle = 67310.54841 \
    + (876600 * 60 * 60 + 8640184.812866)*tut \
    + 0.093104*(tut**2) \
    + (-6.2e-6)*(tut**3)
  return ((gmst_angle % 86400) * w )

def z_rot_ecef(theta, ecef_x_km, ecef_y_km, ecef_z_km):
  Rz = [[math.cos(-theta), -math.sin(-theta), 0], [math.sin(-theta), math.cos(-theta), 0], [0, 0, 1]]
  det_Rz = determinant(Rz)
  adj_Rz = adjugate(Rz)
  Rz_inverse = [[adj_Rz[i][j] / det_Rz for j in range(3)] for i in range(3)]

  eci_x_km = Rz_inverse[0][0]*ecef_x_km + Rz_inverse[0][1]*ecef_y_km + Rz_inverse[0][2]*ecef_z_km
  eci_y_km = Rz_inverse[1][0]*ecef_x_km + Rz_inverse[1][1]*ecef_y_km + Rz_inverse[1][2]*ecef_z_km
  eci_z_km = Rz_inverse[2][0]*ecef_x_km + Rz_inverse[2][1]*ecef_y_km + Rz_inverse[2][2]*ecef_z_km

  return [eci_x_km, eci_y_km, eci_z_km]

# Function to calculate the determinant of a 3x3 matrix
def determinant(matrix):
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]

    det = a*(e*i-f*h) \
        -b*(d*i-f*g) \
        +c*(d*h-e*g)
    
    return det

# Function to calculate the adjugate of a 3x3 matrix
def adjugate(matrix):
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    
    adj = [
        [e*i-f*h, c*h-b*i, b*f-c*e],
        [f*g-d*i, a*i-c*g, c*d-a*f],
        [d*h-e*g, b*g-a*h, a*e-b*d]
    ]

    return adj

# main function
def ecef_to_eci(y, m, d, hr, min, sec, ecef_x_km, ecef_y_km, ecef_z_km):
  jd_frac = ymdhms_tojd(y, m, d, hr, min, sec)
  t_ut = jd_to_tut(jd_frac)
  gmst_angle_rads = tut_to_GMST_rads(t_ut)
  eci = z_rot_ecef(gmst_angle_rads, ecef_x_km, ecef_y_km, ecef_z_km)
  print(eci[0])
  print(eci[1])
  print(eci[2])

  return eci
  
# initialize script arguments
year = None
month = None
day = None
hour = None
min = None
sec = None
ecef_x_km = None
ecef_y_km = None
ecef_z_km = None

# parse script arguments
if len(sys.argv)==10:
  year = int(sys.argv[1])
  month = int(sys.argv[2])
  day = int(sys.argv[3])
  hour = int(sys.argv[4])
  min = int(sys.argv[5])
  sec = float(sys.argv[6])
  ecef_x_km = float(sys.argv[7])
  ecef_y_km = float(sys.argv[8])
  ecef_z_km = float(sys.argv[9])

else:
  print('Usage: python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km')
  exit()

# write script below this line
if __name__ == '__main__':
  ecef_to_eci(year, month, day, hour, min, sec, ecef_x_km, ecef_y_km, ecef_z_km)
else:
  ecef_to_eci(year, month, day, hour, min, sec, ecef_x_km, ecef_y_km, ecef_z_km)