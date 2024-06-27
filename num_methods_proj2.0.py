#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# file_name.py
"""

@author: claudio

Created: Wed Jan 18 16:36:32 2023

Progam utilizes Euler's Method, Improved Euler's Method, and Runge-Kutta Method to numerically solve for 

-----------


"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


print("\n Euler's Method")


def eulers_method(y_0, h):

    def f(x, y): return 3.0*(x) - 6.0*(y)

    global x
    x = np.arange(3.0, 3.60, h)
    Y = np.zeros(len(x))  # creates an array for y
    Y[0] = y_0  # set the first value of the array to y_0

    for i in range(0, len(x) - 1):
        Y[i + 1] = Y[i] + h*f(x[i], Y[i])

        print("\n", "i : ", i, "     ",  # prints the values of i, h, x, y, and k
              "\n h : ", h, "     ",
              "\n x :", np.round(x[i], 2), "     ",
              "\n y : ", np.round(Y[i], 8), "     ",
              "\n k :", np.round(f(x[i], Y[i]), 8))


euler_m = eulers_method(0.1, 0.08)

print("\n Improved Euler's Method")


def improved_eulers_method(y_0, h):

    def f(x, y): return 3.0*(x) - 6.0*(y)

    Z = np.zeros(len(x))  # creates an array for z
    Y1 = np.zeros(len(x))  # creates an array for y
    k1 = np.zeros(len(x))  # creates an array for k1
    k2 = np.zeros(len(x))  # creates an array for k2

    Y1[0] = y_0  # sets the first value of the array to y_0

    for i in range(0, len(x)-1):
        Z[i] = Y1[i] + (h*f(x[i], Y1[i]))  # calculates z
        Y1[i+1] = Y1[i] + \
            (h * ((f(x[i], Y1[i])+f(x[i+1], Z[i]))/2))  # calculates y

        k1[i] = f(x[i], Y1[i])  # calculates k1
        k2[i] = f(x[i+1], Z[i])  # calculates k2

        print("\n", "i : ", i, "     ",  # prints the values of i, h, x, y, z, k1, & k2
              "\n h : ", h, "     ",
              "\n x : ", np.round(x[i], 2), "     ",
              "\n y : ", np.round(Y1[i], 8), "     ",
              "\n k1 : ", np.round(k1[i], 8), "     ",
              "\n z : ", np.round(Z[i], 8), "     ",
              "\n k2 : ", np.round(k2[i], 8))


improved_euler_m = improved_eulers_method(0.1, 0.08)


def runge_kutta_method(y_0, h):

    def f(x, y): return 3.0*(x) - 6.0*(y)

    x = np.arange(3.0, 3.60, h)

    Z2 = np.zeros(len(x))  # creates an array for z2
    Z3 = np.zeros(len(x))  # creates an array for z3
    Z4 = np.zeros(len(x))  # creates an array for z4
    Y2 = np.zeros(len(x))  # creates an array for y
    K1 = np.zeros(len(x))  # creates an array for k1
    K2 = np.zeros(len(x))  # creates an array for k2
    K3 = np.zeros(len(x))  # creates an array for k3
    K4 = np.zeros(len(x))  # creates an array for k4
    Y2[0] = y_0  # sets the first value of the array to y_0

    for i in range(0, len(x)-1):
        K1[i] = f(x[i], Y2[i])  # calculates k1
        K2[i] = (f((x[i]+(h/2)), (Y2[i] + (h*(K1[i]/2)))))  # calculates k2
        K3[i] = (f((x[i]+(h/2)), (Y2[i] + (h*(K2[i]/2)))))  # calculates k3
        K4[i] = (f((x[i]+(h)), (Y2[i] + (K3[i]*h))))  # calculates k4
        Y2[i + 1] = Y2[i] + \
            ((h/6)*((K1[i]+(2*K2[i])+(2*K3[i])+K4[i])))  # calculates y

        Z2[i] = Y2[i] + (K1[i]*(h/2))  # calculates z2
        Z3[i] = Y2[i] + (K2[i]*(h/2))  # calculates z3
        Z4[i] = Y2[i] + (K3[i]*(h/2))  # calculates z4

        print("\n", "i:", i, "     ",  # prints the values of i, h, x, y, z[2,4], & k[1,4]
              "\n h : ", h, "     ",
              "\n x :", np.round(x[i], 2), "     ",
              "\n y : ", np.round(Y2[i], 8), "     ",
              "\n k1 :", np.round(K1[i], 8), "     ",
              "\n z2 :", np.round(Z2[i], 8), "     ",
              "\n k2 :", np.round(K2[i], 8), "     ",
              "\n z3 :", np.round(Z3[i], 8), "     ",
              "\n k3 :", np.round(K3[i], 8), "     ",
              "\n z4 :", np.round(Z4[i], 8), "     ",
              "\n k4 :", np.round(K4[i], 8), "     ",)


print("\n Runge Kutta Method")

runge_kutta_m = runge_kutta_method(0.1, 0.08)
