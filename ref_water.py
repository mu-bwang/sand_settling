# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 21:37:43 2021

@author: wangbinb
"""
import numpy as np

R = 8.31446261815324 # Gas constant, (J/(mol K))

class water():
        """
        Parameters:
        -------------------------------------
                P - Pressure (Pa)
               Pv - vapor pressure (Pa)
                T - Temperature (degree C)
                z - Depth (m)
                S - salinity (psu)
              rho - density (kg/m^3)
               nu - kinematic viscosity (m^2/s)
               mu - dynamic viscosity (Pa s, or N s/m^2)
                k - thermal conductivity (W/(mK))
               cp - heat capacity (J/(kg K))
            sigma - surface tension (N/m)
                c - sound speed (m/s)
               pH - pH value of water
        """
        def __init__(self, temperature = 22, depth = 0, salinity = 0, pH = 8.0):
            self.T = temperature
            self.z = depth
            self.S = salinity
            self.pH = pH
            self.P = self.pressure()
            self.rho = self.density()

            self.nu = self.viscosity()[0]
            self.mu = self.viscosity()[1]
            self.sigma = self.surface_tension()
            self.k = self.thermal_conductivity()
            self.cp = 3997.4
            self.Pv = self.vapor_pressure()
            self.c = self.sound_speed()



        def density(self):

            if self.T < 40:
                # Convert P to bar
                P_in_bar = self.P * 1.e-5

                # Compute the density at atmospheric pressure
                rho_sw_0 = (
                            999.842594 + 6.793952e-2 * self.T - 9.095290e-3 * self.T**2
                            + 1.001685e-4 * self.T**3 - 1.120083e-6 * self.T**4 + 6.536332e-9 * self.T**5
                            + 8.24493e-1 * self.S - 5.72466e-3 * self.S**(3./2.) + 4.8314e-4 * self.S**2
                            - 4.0899e-3 * self.T*self.S + 7.6438e-5 * self.T**2 * self.S - 8.2467e-7 * self.T**3 *
                            self.S + 5.3875e-9 * self.T**4 * self.S + 1.0227e-4 * self.T * self.S**(3./2.)
                            - 1.6546e-6 * self.T**2 * self.S**(3./2.)
                            )

                # Compute the pressure correction coefficient
                K = (
                     19652.21 + 148.4206 * self.T - 2.327105 * self.T**2 + 1.360477e-2 * self.T**3
                     - 5.155288e-5 * self.T**4 + 3.239908 * P_in_bar + 1.43713e-3 * self.T * P_in_bar
                     + 1.16092e-4 * self.T**2 * P_in_bar - 5.77905e-7 * self.T**3 * P_in_bar
                     + 8.50935e-5 * P_in_bar**2 - 6.12293e-6 * self.T * P_in_bar**2
                     + 5.2787e-8 * self.T**2 * P_in_bar**2 + 54.6746 * self.S - 0.603459 * self.T * self.S
                     + 1.09987e-2 * self.T**2 * self.S - 6.1670e-5 * self.T**3 * self.S
                     + 7.944e-2 * self.S**(3./2.) + 1.64833e-2 * self.T * self.S**(3./2.)
                     - 5.3009e-4 * self.T**2 * self.S**(3./2.) + 2.2838e-3 * P_in_bar * self.S
                     - 1.0981e-5 * self.T * P_in_bar * self.S - 1.6078e-6 * self.T**2 * P_in_bar * self.S
                     + 1.91075e-4 * P_in_bar * self.S**(3./2.) - 9.9348e-7 * P_in_bar**2 * self.S
                     + 2.0816e-8 * self.T * P_in_bar**2 * self.S + 9.1697e-10 * self.T**2 * P_in_bar**2 * self.S
                     )

                rho = rho_sw_0 / (1 - P_in_bar / K)

            else:
                # Convert P to MPa
                P_in_MPa = self.P / 1.e6

                # Summations
                left_col = 9.9920571e2 + 9.5390097e-2 * self.T - 7.6186636e-3 * self.T**2 + \
                    3.1305828e-5 * self.T**3 - 6.1737704e-8  * self.T**4 + 4.3368858e-1 * P_in_MPa + \
                    2.5495667e-5 * P_in_MPa*self.T**2 - 2.8988021e-7 * P_in_MPa*self.T**3 + \
                    9.5784313e-10 * P_in_MPa*self.T**4 + 1.7627497e-3 * P_in_MPa**2 - 1.2312703e-4 * P_in_MPa**2*self.T \
                    + 1.3659381e-6 * P_in_MPa**2*self.T**2 + 4.0454583e-9 * P_in_MPa**2*self.T**3 - 1.4673241e-5 \
                    * P_in_MPa**3 + 8.8391585e-7 * P_in_MPa**3*self.T - 1.1021321e-9 * P_in_MPa**3*self.T**2 + \
                    4.2472611e-11 * P_in_MPa**3*self.T**3 - 3.9591772e-14 * P_in_MPa**3*self.T**4

                right_col = -7.99992230e-1 * self.S + 2.40936500e-3 * self.S*self.T - 2.58052775e-5 * \
                    self.S*self.T**2 + 6.85608405e-8 * self.S*self.T**3 + 6.29761106e-4 * P_in_MPa*self.S - \
                    9.36263713e-7 * P_in_MPa**2*self.S

                rho = left_col - right_col
            return rho

        def viscosity(self):
            # mu = 2.414e-5 * 10 ** (247.8 / (self.T+273.15-140))


            # Get the fit coefficients
            a = np.array([1.5700386464E-01, 6.4992620050E+01, -9.1296496657E+01,
                          4.2844324477E-05, 1.5409136040E+00, 1.9981117208E-02,
                          -9.5203865864E-05, 7.9739318223E+00, -7.5614568881E-02,
                          4.7237011074E-04])

            # Compute the viscosity of pure water at given temperature
            mu_w = a[3] + 1./(a[0] * (self.T + a[1])**2 + a[2])

            # Correct for salinity
            S = self.S / 1000.
            A = a[4] + a[5] * self.T + a[6] * self.T**2
            B = a[7] + a[8] * self.T + a[9] * self.T**2
            mu_0 = mu_w * (1. + A * S + B * S**2)

            # And finally for pressure
            P = self.P * 0.00014503773800721815
            mu = mu_0 * (0.9994 + 4.0295e-5 * P + 3.1062e-9 * P**2)
            nu = mu/ self.rho

            return nu, mu

        def pressure(self):
            P = 1.013e5 * (1 + 0.1 * self.z)
            return P

        def surface_tension(self):
            # Equations in Sharqawy using deg C and g/kg
            S = self.S / 1000.

            # Use equations (27) for pure water surface tension (N/m)
            sigma_w = 0.2358 * (1. - (self.T + 273.15) / 647.096)**1.256 * (1. - 0.625 *
                      (1. - (self.T + 273.15) / 647.096))

            # Equation (28) gives the salinity correction
            if self.T < 40:
                # Salinity correction only valid [0, 40] deg C
                sigma = sigma_w * (1. + (0.000226 * self.T + 0.00946) * np.log(1. + 0.0331 *
                        S))
            else:
                # No available salinity correction for hot cases
                sigma = sigma_w

            return sigma

        def thermal_conductivity(self):
            # Thermal conductivity equations use T_68 in deg C
            T_68 = ((self.T+273.15) - 0.0682875) / (1.0 - 0.00025)
            T_68  = T_68 - 273.15

            # Salinity is in g/kg and pressure in MPa
            S = self.S / 1000.
            P = self.P * 1e-6

            # Compute the thermal conductivity from Table 4
            if T_68 < 30.:
                # Equation (15)
                k = 0.55286 + 3.4025e-4 * P + 1.8364e-3 * T_68 - 3.3058e-7 * \
                       T_68**3

            else:
                # Equation (13)
                k = 10.**(np.log10(240. + 0.0002 * S) + 0.434 * (2.3 - (343.5 +
                       0.037 * S) / (T_68 + 273.15)) * (1. - (T_68 + 273.15) /
                       (647. + 0.03 * S)) ** 0.333) / 1000.

            return k

        def vapor_pressure(self):
            Pv = (0.61094 * np.exp(17.625*self.T/(self.T+243.04)))*1e3;
            return Pv

        def sound_speed(self):
            c = 1449.2 + 4.6 * self.T - 0.055 * self.T**2 + 0.00029 * self.T**3 \
                + (1.34 - 0.01 * self.T) * (self.S - 35) + 0.016 * self.z
            return c
