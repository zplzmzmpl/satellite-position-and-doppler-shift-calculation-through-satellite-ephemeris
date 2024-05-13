import numpy as np
import math
from math import sqrt

def cal_n(sqrta,delta_n):
  n0 = sqrt(3.986005e14) / np.power(sqrta,3)
  return n0 + delta_n

def cal_M(M0, t, te, n):
  return M0 + n*(t-te)

def cal_E(M,e):
  count = 0
  E0 = 0
  E = 1
  while abs(E-E0)>1e-14:
    count += 1
    E = E0
    E0 = M + e * np.sin(E)
    if count > 3000:
      print('error: too many iterations!')
      break
  return E

def cal_f(E,e):
  cosf = (np.cos(E) - e) / (1 - e * np.cos(E))
  sinf = sqrt(1 - np.power(e,2)) * np.sin(E) / (1 - e * np.cos(E))
  return np.arctan2(sinf,cosf)

def cal_uk(f,w0):
  return f + w0

def cal_sigma_u_r_i(Cuc, Cus, Crc, Crs, Cic, Cis, uk):
  return Cuc*np.cos(2*uk) + Cus*np.sin(2*uk), Crc*np.cos(2*uk) + Crs*np.sin(2*uk), Cic*np.cos(2*uk) + Cis*np.sin(2*uk)

def cal_uk_rk_ik(uk, sqrta, e, E, sigma_u, sigma_r, sigma_i, i0, t, te, delta_d):
  return uk + sigma_u, np.power(sqrta,2)*(1-e*np.cos(E)) + sigma_r, i0 + sigma_i + delta_d*(t-te)

def cal_x_y(uk, rk):
  return rk*np.cos(uk), rk*np.sin(uk)

def cal_longitude_for_GEO(omega0, omegaDot, t, te, we=7.2921151467e-5):
  return omega0 + omegaDot*(t-te) - we*te

def cal_longitude_for_else(omega0, omegaDot, t, te, we=7.2921151467e-5):
  return omega0 + omegaDot*(t-te) - we*t

def cal_X_Y_Z(x0,y0,L,i):
  x = x0*np.cos(L) - y0*np.cos(i)*np.sin(L)
  y = x0*np.sin(L) + y0*np.cos(i)*np.cos(L)
  z = y0*np.sin(i)
  return x,y,z

def cal_GEO_coordination(t, te, x, y, z, we=7.2921151467e-5):
  phi = - (5 * np.pi / 180)
  theta = we * (t - te)
  Rx_phi = np.array([[1,0,0],
                    [0,np.cos(phi),np.sin(phi)],
                    [0,-np.sin(phi),np.cos(phi)]])
  Rz = np.array([[np.cos(theta),np.sin(theta),0],
                 [-np.sin(theta),np.cos(theta),0],
                 [0,0,1]])
  return Rz.dot(Rx_phi).dot(np.array([x,y,z]))

def run_cal_coordinations(path, t):
  with open(path, 'r') as file:
      # read by line
      lines = file.readlines()

  # create a dict
  variables = {}

  for line in lines:
      # split key value pair
      key, value = line.strip().split(': ')
      # add value
      variables[key] = float(value)

  id = variables['svid']
  te = variables['toe']
  e = variables['e']
  sqrta = variables['sqrta']
  M0 = variables['M0']
  w0 = variables['w0']
  i0 = variables['i0']
  omega0 = variables['Omega0']
  delta_d = variables['di/dt']
  delta_n = variables['deltaN']
  omegaDot = variables['OmegaDot']
  Cuc = variables['Cuc']
  Cus = variables['Cus']
  Crc = variables['Crc']
  Crs = variables['Crs']
  Cic = variables['Cic']
  Cis = variables['Cis']

  del variables, lines

  n = cal_n(sqrta,delta_n)
  M = cal_M(M0, t, te, n)
  E = cal_E(M,e)
  f = cal_f(E,e)
  uk = cal_uk(f,w0)
  sigma_u, sigma_r, sigma_i = cal_sigma_u_r_i(Cuc, Cus, Crc, Crs, Cic, Cis, uk)
  uk, rk, ik = cal_uk_rk_ik(uk, sqrta, e, E, sigma_u, sigma_r, sigma_i, i0, t, te, delta_d)
  x0, y0 = cal_x_y(uk, rk)

  if 101 <= id <= 105 or id == 159:
    L = cal_longitude_for_GEO(omega0,omegaDot,t,te)
    X, Y, Z = cal_X_Y_Z(x0, y0, L, ik)
    # print('calculated X,Y,Z:',X,Y,Z)
    x, y, z = cal_GEO_coordination(t, te, X, Y, Z)
    print('calculated x,y,z for GEO satellite:',x,y,z)
    return x, y, z
  else:
    L = cal_longitude_for_else(omega0,omegaDot,t,te)
    X, Y, Z = cal_X_Y_Z(x0, y0, L, ik)
    print('calculated X,Y,Z:',X,Y,Z)
    return X, Y, Z

def calculate_omega_and_vs(T, orbit_radius):
  omega = 2*np.pi/T
  vs = orbit_radius*1000*omega
  return omega, vs

def calculate_doppler_shift(x1, y1, z1, x2, y2, z2, xr, yr, zr, L1 = 1575.42, c = 3e8):
  # Calculate d1 and d2
  d1 = sqrt((x1 - xr)**2 + (y1 - yr)**2 + (z1 - zr)**2)
  d2 = sqrt((x2 - xr)**2 + (y2 - yr)**2 + (z2 - zr)**2)
  print('d1:',d1,'d2:',d2)
  # Calculate vd
  vd = abs(d1 - d2)
  # Calculate fd (Doppler frequency shift)
  fd = vd * L1 / c
  return fd

def run_doppler(path, time, xr, yr, zr):
  for t in range(time,time+3):
    x1, y1, z1 = run_cal_coordinations(path, t)
    x2, y2, z2 = run_cal_coordinations(path, t+1)
    fd = calculate_doppler_shift(x1, y1, z1, x2, y2, z2, xr, yr, zr)
    print("Doppler frequency shift:", fd)


  

if __name__ == '__main__':

  # part 1
  path = r'C:\PlatformData\beidou\GNSS_DATA\Eph2.txt'
  t = 439075
  run_cal_coordinations(path, t)

  # part 2
  # GEO_T = 23*3600 + 56*60 + 4
  # GEO_radius = 42160 # kilometer
  # GEO_omega, GEO_vs = calculate_omega_and_vs(GEO_T, GEO_radius)
  # print('omeg and average tangential velocity of GEO satellite:', GEO_omega, GEO_vs)

  # IGSO_T = 23*3600 + 56*60 + 4
  # IGSO_radius = 42270
  # IGSO_omega, IGSO_vs = calculate_omega_and_vs(IGSO_T, IGSO_radius)
  # print('omeg and average tangential velocity of IGSO satellite:', IGSO_omega, IGSO_vs)

  # MEO_T = 11*3600 + 58*60 + 2
  # MEO_radius = 27980
  # MEO_omega, MEO_vs = calculate_omega_and_vs(MEO_T, MEO_radius)
  # print('omeg and average tangential velocity of MEO satellite:', MEO_omega, MEO_vs)

  # x1, y1, z1 = -14440511.39, 21944026.56, -2229631.20 # GPS 02 MEO	439074
  # x2, y2, z2 = -14440946.36, 21944131.10, -2226438.89 # GPS 02 MEO	439075
  # x2, y2, z2 = -14441381.09, 21944235.34, -2223246.53 # GPS 02 MEO	439076
  # x2, y2, z2 = -14440810.83, 21945321.20, -2217421.51 # GPS 02 MEO	482156
  # xr, yr, zr = -2291885.56, 5002448.07, 3214696.11
  # run_doppler(path, 439074, xr, yr, zr)

  # x1, y1, z1 = -14461105.87, 21948708.31, -2076346.88 # GPS 02 MEO	439122
  # x2, y2, z2 = -14461528.77, 21948798.53, -2073152.32 # GPS 02 MEO	439123
  # x2, y2, z2 = -14461951.41, 21948888.46, -2069957.72 # GPS 02 MEO	439124
  # x2, y2, z2 = -14461360.43, 21949948.30, -2064129.98 # GPS 02 MEO	482204
  xr, yr, zr = -2291885.26, 5002448.23, 3214695.95 # Receiver coordinates
  run_doppler(path, 439122, xr, yr, zr)

  # x1, y1, z1 = -14148225.52, 38383202.86, 9658585.11 # BDS 106 IGSO 439074
  # x2, y2, z2 = -14147357.00, 38384120.02, 9656182.01 # BDS 106 IGSO 439075
  # x2, y2, z2 = -14146488.35, 38385037.04, 9653778.85 # BDS 106 IGSO 439076
  # x2, y2, z2 = -14104116.40, 38387940.78, 9704742.90 # BDS 106 IGSO 525238
  # xr, yr, zr = -2291885.56, 5002448.07, 3214696.11 # Receiver coordinates
  # run_doppler(path, 439074, xr, yr, zr)

  # x1, y1, z1 = -34318296.01, 24481025.59,-1113665.51 # BDS 101 GEO 439074
  # x2, y2, z2 = -34318297.57, 24481026.88, -1113616.78 # BDS 101 GEO 439075
  # x2, y2, z2 = -34318299.14, 24481028.16, -1113568.04 # BDS 101 GEO 439076
  # x2, y2, z2 = -34324482.43, 24472246.56, -1116047.11 # BDS 101 GEO 525238
  # xr, yr, zr = -2291885.56, 5002448.07, 3214696.11 # Receiver coordinates
  # run_doppler(path, 439074, xr, yr, zr)

  # t = 1 # Time interval
  # L1 = 1575.42 # Wave frequency MHz
  # c = 3e8 # Speed of light in m/s

  
  