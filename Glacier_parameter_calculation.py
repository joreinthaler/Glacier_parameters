# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 12:19:14 2023

@author: jreintha
"""
#The following script calculates glacier parameters for different time points (depending on the users inventory).
#In this case, glacier parameters from outlines from the Little Ice Age (time around 1850), the Randolph glacier inventory (time around 2003), and 2015 (here called cci). Furthermore other datasets from the 60s and 70s were included (here called SxS)
#The script follows the methodology of the follofing references: 
#   Haeberli & Hoelzle 1995 (10.3189/s0260305500015834) 
#   Hoelzle et al. 2003 (10.1016/S0921-8181(02)00223-0) 
#   Hoelzle et al. 2007 (10.1016/j.gloplacha.2006.07.001)

#The script also compares the results to volume changes from an external GIS-based reconstruction. This can be removed if only the glacier parameters want to be generated.

import os
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from matplotlib.lines import Line2D
import seaborn as sns

os.chdir('filepath')   

#%% Define constants
f = 0.8  # shape factor
p = 900  # density of ice in kg/m^3
g = 9.81  # gravitational acceleration in m/s^2
n = 3  # a realistic value
A = 0.16  # in a^1 bar^3

#%% Import Excel file and create a pandas DataFrame object
file_path_LIA = "filepath"
LIA_df = pd.read_excel(file_path_LIA)
file_path_cci = "filepath.xlsx"
cci_df = pd.read_excel(file_path_cci)
file_path_RGI = "filepath.xlsx"
RGI_df = pd.read_excel(file_path_RGI)
file_path_SxS = "filepath.xlsx"
SxS_df = pd.read_excel(file_path_SxS)

file_path_Flow_velocity = "filepath.xlsx"
flow_df = pd.read_excel(file_path_Flow_velocity)

xls = pd.ExcelFile("filepath.xlsx")
Parameter_df = xls.parse("Sheet1")
GIS_df = xls.parse("Sheet2")
GIS_region = xls.parse("Sheet3")
Missing_gl = pd.ExcelFile("filepath.xlsx")
Missing_gl_df = Missing_gl.parse("V_HaHo")

#  Part I: Calculation of glacier parameters
# %% Calcualte variables LIA
# Calculate Hmean and dH variables
LIA_df["Hmean"] = (LIA_df["MAX"] + LIA_df["MIN"]) / 2
LIA_df["dH"] = LIA_df["MAX"] - LIA_df["MIN"]

# LIA_df['dH_km'] = (LIA_df['MAX'] - LIA_df['MIN'])/1000
LIA_df["dHa"] = (LIA_df["MAX"] - LIA_df["MIN"]) / 2

# area in km2
LIA_df["AREA"] = LIA_df["AREA"] * 0.000001

# Calculate the length
LIA_df["L0"] = LIA_df["LENGTH"]
LIA_df["L0_m"] = LIA_df["L0"] * 1000

# Length of ablation area (La)
LIA_df["La"] = np.where(LIA_df["L0"] < 2, LIA_df["L0"] * 0.5, LIA_df["L0"] * 0.75)
LIA_df["La_m"] = LIA_df["La"] * 1000

# Calculate average slope

LIA_df["a"] = LIA_df.apply(
    lambda row: math.degrees(math.atan(row["dH"] / row["L0_m"])), axis=1
)

##LIA_df['dH'] / LIA_df['L0'].apply(lambda x: math.sqrt(1 + x**2))

# Calculate average slope of ablation area
LIA_df["aa"] = LIA_df.apply(
    lambda row: math.degrees(math.atan(row["dHa"] / row["La_m"])), axis=1
)

# Calculate the basal shear stress
# threshold
LIA_df["tf"] = (
    0.005 + 1.598 * (LIA_df["dH"] / 1000) - 0.435 * ((LIA_df["dH"] / 1000) ** 2)
)
LIA_df.loc[LIA_df["dH"] > 1600, "tf"] = 1.5
LIA_df["tf_Pa"] = LIA_df["tf"] * 10 ** 5

# Calculate sinus alpha
def sine(x):
    return math.sin(math.radians(x))


LIA_df["sin_a"] = LIA_df["a"].apply(sine)
LIA_df["sin_aa"] = LIA_df["aa"].apply(sine)

# Calculate average ice depth along flowline (with or without f??)
LIA_df["hf"] = LIA_df["tf_Pa"] / (f * p * g * LIA_df["sin_a"])
LIA_df["hfa"] = LIA_df["tf_Pa"] / (f * p * g * LIA_df["sin_aa"])

# Calculate average thickness
LIA_df["hF"] = (math.pi / 4) * LIA_df["hf"]

# Calculate maximum thickness
LIA_df["h_max"] = 2.5 * LIA_df["hfa"]

# Calculate Volume
LIA_df["V"] = LIA_df["AREA"] * 1000000 * LIA_df["hF"]
LIA_df["V_km3"] = LIA_df["V"] / 1e9

# ------------------------------------------------------------------------------
# Mass Balanca
# Calclate mass balance gradient (0.75 m per 100m) in ablation area
LIA_df["db/dH"] = -0.75 / 100

# Calculate mass balance at glacier thongue
LIA_df["bt"] = (LIA_df["db/dH"] * (LIA_df["Hmean"] - LIA_df["MIN"])) * -1

# ------------------------------------------------------------------------------
# Flow velocity
# Flow velocity along the flowline in the lower half of the ablation area
LIA_df["u_ma"] = ((3 * LIA_df["bt"] / 4) * (LIA_df["La_m"] / 2)) / LIA_df["hfa"]

# Surface veolociy
LIA_df["u_sa"] = LIA_df["u_ma"]

# Calculate component of ice deformation
LIA_df["u_da"] = 2 * A * (LIA_df["tf"] ** n) * LIA_df["hfa"] / (n + 1)

# Sliding velocity in the ablation area
LIA_df["u_ba"] = LIA_df["u_sa"] - LIA_df["u_da"]

# ------------------------------------------------------------------------------
# Responce time
# Calculate responce time
LIA_df["t_resp"] = LIA_df["h_max"] / LIA_df["bt"]

# Calculate kinematic wave velocity
LIA_df["c"] = 4 * LIA_df["u_sa"]

# Calculate reaction time
LIA_df["t_react"] = LIA_df["La_m"] / LIA_df["c"]

#Calculate glacier volume
LIA_df['V_km3'].replace([np.inf, -np.inf, np.nan], 0, inplace=True)
LIA_volume = LIA_df['V_km3'].sum()
print('LIA ice volume'), print(LIA_volume)

#Calculate glacier volume per region
LIA_sum_vol_region = LIA_df.groupby('region_id')['V_km3'].sum()
print("Sum of 'V_km3' per region in LIA_df:")
print(LIA_sum_vol_region)

#Calculate sum of area
LIA_df["AREA"].replace([np.inf, -np.inf, np.nan], 0, inplace=True)
LIA_area = LIA_df["AREA"].sum()

# Calculate regional sums and averages
region_summary_LIA = (
    LIA_df.groupby("region_id")
    .agg({"V_km3": "sum", "AREA": "sum", "MIN": "mean"})
    .reset_index()
)

# Rename the columns for clarity
region_summary_LIA.rename(
    columns={
        "V_km3_LIA": "Total_V_km3_LIA",
        "AREA_LIA": "Total_AREA_LIA",
        "MIN_LIA": "Average_MIN_LIA",
    },
    inplace=True,
)

# %% Calcualte variables cci 2015

# Calculate Hmean and dH variables
cci_df["Hmean"] = (cci_df["MAX"] + cci_df["MIN"]) / 2
cci_df["dH"] = cci_df["MAX"] - cci_df["MIN"]

# cci_df['dH_km'] = (cci_df['MAX'] - cci_df['MIN'])/1000
cci_df["dHa"] = (cci_df["MAX"] - cci_df["MIN"]) / 2

# area in km2
cci_df["AREA"] = cci_df["AREA"] * 0.000001

# Calculate the length (km)
cci_df["L0"] = cci_df["LENGTH"]
# Calculate the length (m)
cci_df["L0_m"] = cci_df["L0"] * 1000

# Length of ablation area (La)
cci_df["La"] = np.where(cci_df["L0"] < 2, cci_df["L0"] * 0.5, cci_df["L0"] * 0.75)
cci_df["La_m"] = cci_df["La"] * 1000
# Calculate average slope

cci_df["a"] = cci_df.apply(
    lambda row: math.degrees(math.atan(row["dH"] / row["L0_m"])), axis=1
)

##cci_df['dH'] / cci_df['L0'].apply(lambda x: math.sqrt(1 + x**2))

# Calculate average slope of ablation area
cci_df["aa"] = cci_df.apply(
    lambda row: math.degrees(math.atan(row["dHa"] / row["La_m"])), axis=1
)

# Calculate the basal shear stress
# threshold
cci_df["tf"] = (
    0.005 + 1.598 * (cci_df["dH"] / 1000) - 0.435 * ((cci_df["dH"] / 1000) ** 2)
)
cci_df.loc[cci_df["dH"] > 1600, "tf"] = 1.5
cci_df["tf_Pa"] = cci_df["tf"] * 10 ** 5

# Calculate sinus alpha


def sine(x):
    return math.sin(math.radians(x))


cci_df["sin_a"] = cci_df["a"].apply(sine)
cci_df["sin_aa"] = cci_df["aa"].apply(sine)


# Calculate average ice depth along flowline (with or without f??)
cci_df["hf"] = cci_df["tf_Pa"] / (f * p * g * cci_df["sin_a"])
cci_df["hfa"] = cci_df["tf_Pa"] / (f * p * g * cci_df["sin_aa"])

# Calculate average thickness
cci_df["hF"] = (math.pi / 4) * cci_df["hf"]

# Calculate maximum thickness
cci_df["h_max"] = 2.5 * cci_df["hfa"]

# Calculate Volume
cci_df["V"] = cci_df["AREA"] * 1000000 * cci_df["hF"]
cci_df["V_km3"] = cci_df["V"] / 1e9

# ------------------------------------------------------------------------------
# Mass Balanca
# Calclate mass balance gradient (0.75 m w.e. per 100m) in ablation area
cci_df["db/dH"] = 0.75 / 100

# Calculate mass balance at glacier thongue
cci_df["bt"] = cci_df["db/dH"] * (cci_df["Hmean"] - cci_df["MIN"])

# ------------------------------------------------------------------------------
# Flow velocity
# Flow velocity along the flowline in the lower half of the ablation area
cci_df["u_ma"] = ((3 * cci_df["bt"] / 4) * (cci_df["La_m"] / 2)) / cci_df["hfa"]

# Surface veolociy
cci_df["u_sa"] = cci_df["u_ma"]

# Calculate component of ice deformation
cci_df["u_da"] = 2 * A * (cci_df["tf"] ** n) * cci_df["hfa"] / (n + 1)

# Sliding velocity in the ablation area
cci_df["u_ba"] = cci_df["u_sa"] - cci_df["u_da"]

# ------------------------------------------------------------------------------
# Responce time
# Calculate responce time
cci_df["t_resp"] = cci_df["h_max"] / cci_df["bt"]

# Calculate kinematic wave velocity
cci_df["c"] = 4 * cci_df["u_sa"]

# Calculate reaction time
cci_df["t_react"] = cci_df["La_m"] / cci_df["c"]

#Calculate glacier volume
cci_volume = cci_df.loc[cci_df['LIA_ID'] < 5000, 'V_km3'].sum()
cci_volume_all = cci_df['V_km3'].sum()
print('cci ice volume'), print(cci_volume)
print('cci ice volume_all'), print(cci_volume_all)

#Calculate glacier volume per region
cci_sum_vol_region = cci_df.groupby('region_id')['V_km3'].sum()
print("\nSum of 'V_km3' per region in cci_df:")
print(cci_sum_vol_region)

#Calculate sum of area
cci_area = cci_df.loc[cci_df["LIA_ID"] < 5000, "AREA"].sum()
cci_area_all = cci_df["AREA"].sum()

#Coverage
coverage_cci = 100 / cci_area_all * cci_area
print("Coverage cci"), print(coverage_cci)

# Calculate regional sums and averages
region_summary_cci = (
    cci_df.groupby("region_id")
    .agg({"V_km3": "sum", "AREA": "sum", "MIN": "mean"})
    .reset_index()
)

# Rename the columns for clarity
region_summary_cci.rename(
    columns={
        "V_km3_LIA": "Total_V_km3_LIA",
        "AREA_LIA": "Total_AREA_LIA",
        "MIN_LIA": "Average_MIN_LIA",
    },
    inplace=True,
)

# %% Calcualte variables RGI 2003
# Calculate Hmean and dH variables
RGI_df["Hmean"] = (RGI_df["MAX"] + RGI_df["MIN"]) / 2
RGI_df["dH"] = RGI_df["MAX"] - RGI_df["MIN"]

# RGI_df['dH_km'] = (RGI_df['MAX'] - RGI_df['MIN'])/1000
RGI_df["dHa"] = (RGI_df["MAX"] - RGI_df["MIN"]) / 2

# area in km2
RGI_df["AREA"] = RGI_df["AREA"] * 0.000001

# Calculate the length
RGI_df["L0_m"] = RGI_df["LENGTH"]
RGI_df["L0"] = RGI_df["L0_m"] / 1000

# Length of ablation area (La)
RGI_df["La"] = np.where(RGI_df["L0"] < 2, RGI_df["L0"] * 0.5, RGI_df["L0"] * 0.75)
RGI_df["La_m"] = RGI_df["La"] * 1000
# Calculate average slope

RGI_df["a"] = RGI_df.apply(
    lambda row: math.degrees(math.atan(row["dH"] / row["L0_m"])), axis=1
)

##RGI_df['dH'] / RGI_df['L0'].apply(lambda x: math.sqrt(1 + x**2))

# Calculate average slope of ablation area
RGI_df["aa"] = RGI_df.apply(
    lambda row: math.degrees(math.atan(row["dHa"] / row["La_m"])), axis=1
)

# Calculate the basal shear stress
# threshold
RGI_df["tf"] = (
    0.005 + 1.598 * (RGI_df["dH"] / 1000) - 0.435 * ((RGI_df["dH"] / 1000) ** 2)
)
RGI_df.loc[RGI_df["dH"] > 1600, "tf"] = 1.5
RGI_df["tf_Pa"] = RGI_df["tf"] * 10 ** 5

# Calculate sinus alpha
def sine(x):
    return math.sin(math.radians(x))


RGI_df["sin_a"] = RGI_df["a"].apply(sine)
RGI_df["sin_aa"] = RGI_df["aa"].apply(sine)

# Calculate average ice depth along flowline (with or without f??)
RGI_df["hf"] = RGI_df["tf_Pa"] / (f * p * g * RGI_df["sin_a"])
RGI_df["hfa"] = RGI_df["tf_Pa"] / (f * p * g * RGI_df["sin_aa"])

# Calculate average thickness
RGI_df["hF"] = (math.pi / 4) * RGI_df["hf"]

# Calculate maximum thickness
RGI_df["h_max"] = 2.5 * RGI_df["hfa"]

# Calculate Volume
RGI_df["V"] = RGI_df["AREA"] * 1000000 * RGI_df["hF"]
RGI_df["V_km3"] = RGI_df["V"] / 1e9
# Calculate the basal shear stress

# ------------------------------------------------------------------------------
# Mass Balanca
# Calclate mass balance gradient (0.75 m per 100m) in ablation area
RGI_df["db/dH"] = 0.75 / 100

# Calculate mass balance at glacier thongue
RGI_df["bt"] = RGI_df["db/dH"] * (RGI_df["Hmean"] - RGI_df["MIN"])

# ------------------------------------------------------------------------------
# Flow velocity
# Flow velocity along the flowline in the lower half of the ablation area
RGI_df["u_ma"] = ((3 * RGI_df["bt"] / 4) * (RGI_df["La_m"] / 2)) / RGI_df["hfa"]

# Surface veolociy
RGI_df["u_sa"] = RGI_df["u_ma"]

# Calculate component of ice deformation
RGI_df["u_da"] = 2 * A * (RGI_df["tf"] ** n) * RGI_df["hfa"] / (n + 1)

# Sliding velocity in the ablation area
RGI_df["u_ba"] = RGI_df["u_sa"] - RGI_df["u_da"]

# ------------------------------------------------------------------------------
# Responce time
# Calculate responce time
RGI_df["t_resp"] = RGI_df["h_max"] / RGI_df["bt"]

# Calculate kinematic wave velocity
RGI_df["c"] = 4 * RGI_df["u_sa"]

# Calculate reaction time
RGI_df["t_react"] = RGI_df["La_m"] / RGI_df["c"]

#Calculate glacier volume
RGI_volume = RGI_df.loc[RGI_df['LIA_ID'] < 5000, 'V_km3'].sum()
RGI_volume_all = RGI_df['V_km3'].sum()
print('RGI ice volume'), print(RGI_volume)
print('RGI ice volume all'), print(RGI_volume_all)

#Calculate glacier volume per region
RGI_sum_vol_region = RGI_df.groupby('region_id')['V_km3'].sum()
print("\nSum of 'V_km3' per region in cci_df:")
print(RGI_sum_vol_region)

#Calculate sum of area
RGI_area = RGI_df.loc[RGI_df["LIA_ID"] < 5000, "AREA"].sum()
RGI_area_all = RGI_df["AREA"].sum()

# Coverage
coverage_RGI = 100 / RGI_area_all * RGI_area
print("Coverage RGI"), print(coverage_RGI)

# Calculate regional sums and averages
region_summary_RGI = (
    cci_df.groupby("region_id")
    .agg({"V_km3": "sum", "AREA": "sum", "MIN": "mean"})
    .reset_index()
)

# Rename the columns for clarity
region_summary_RGI.rename(
    columns={
        "V_km3_LIA": "Total_V_km3_LIA",
        "AREA_LIA": "Total_AREA_LIA",
        "MIN_LIA": "Average_MIN_LIA",
    },
    inplace=True,
)
# %% Calcualte variables SxS 1960s and 70s
# Calculate Hmean and dH variables
SxS_df["Hmean"] = (SxS_df["MAX"] + SxS_df["MIN"]) / 2
SxS_df["dH"] = SxS_df["MAX"] - SxS_df["MIN"]

# SxS_df['dH_km'] = (SxS_df['MAX'] - SxS_df['MIN'])/1000
SxS_df["dHa"] = (SxS_df["MAX"] - SxS_df["MIN"]) / 2

# area in km2
SxS_df["AREA"] = SxS_df["AREA"] * 0.000001

# Calculate the length
SxS_df["L0"] = SxS_df["LENGTH"]
SxS_df["L0_m"] = SxS_df["L0"] * 1000

# Length of ablation area (La)
SxS_df["La"] = np.where(SxS_df["L0"] < 2, SxS_df["L0"] * 0.5, SxS_df["L0"] * 0.75)
SxS_df["La_m"] = SxS_df["La"] * 1000
# Calculate average slope

SxS_df["a"] = SxS_df.apply(
    lambda row: math.degrees(math.atan(row["dH"] / row["L0_m"])), axis=1
)

##SxS_df['dH'] / SxS_df['L0'].apply(lambda x: math.sqrt(1 + x**2))

# Calculate average slope of ablation area
SxS_df["aa"] = SxS_df.apply(
    lambda row: math.degrees(math.atan(row["dHa"] / row["La_m"])), axis=1
)

# Calculate the basal shear stress
# threshold
SxS_df["tf"] = (
    0.005 + 1.598 * (SxS_df["dH"] / 1000) - 0.435 * ((SxS_df["dH"] / 1000) ** 2)
)
SxS_df.loc[SxS_df["dH"] > 1600, "tf"] = 1.5
SxS_df["tf_Pa"] = SxS_df["tf"] * 10 ** 5

# Calculate sinus alpha


def sine(x):
    return math.sin(math.radians(x))


SxS_df["sin_a"] = SxS_df["a"].apply(sine)
SxS_df["sin_aa"] = SxS_df["aa"].apply(sine)


# Calculate average ice depth along flowline (with or without f??)
SxS_df["hf"] = SxS_df["tf_Pa"] / (f * p * g * SxS_df["sin_a"])
SxS_df["hfa"] = SxS_df["tf_Pa"] / (f * p * g * SxS_df["sin_aa"])

# Calculate average thickness
SxS_df["hF"] = (math.pi / 4) * SxS_df["hf"]

# Calculate maximum thickness
SxS_df["h_max"] = 2.5 * SxS_df["hfa"]

# Calculate Volume
SxS_df["V"] = SxS_df["AREA"] * 1000000 * SxS_df["hF"]
SxS_df["V_km3"] = SxS_df["V"] / 1e9

# ------------------------------------------------------------------------------
# Mass Balance
# Calclate mass balance gradient (0.75 m per 100m) in ablation area
SxS_df["db/dH"] = 0.75 / 100

# Calculate mass balance at glacier thongue
SxS_df["bt"] = SxS_df["db/dH"] * (SxS_df["Hmean"] - SxS_df["MIN"])

# ------------------------------------------------------------------------------
# Flow velocity
# Flow velocity along the flowline in the lower half of the ablation area
SxS_df["u_ma"] = ((3 * SxS_df["bt"] / 4) * (SxS_df["La_m"] / 2)) / SxS_df["hfa"]

# Surface veolociy
SxS_df["u_sa"] = SxS_df["u_ma"]

# Calculate component of ice deformation
SxS_df["u_da"] = 2 * A * (SxS_df["tf"] ** n) * SxS_df["hfa"] / (n + 1)

# Sliding velocity in the ablation area
SxS_df["u_ba"] = SxS_df["u_sa"] - SxS_df["u_da"]

# ------------------------------------------------------------------------------
# Responce time
# Calculate responce time
SxS_df["t_resp"] = SxS_df["h_max"] / SxS_df["bt"]

# Calculate kinematic wave velocity
SxS_df["c"] = 4 * SxS_df["u_sa"]

# Calculate reaction time
SxS_df["t_react"] = SxS_df["La_m"] / SxS_df["c"]


# %% Calcualte variables Missing glaciers
# Calculate Hmean and dH variables
Missing_gl_df["Hmean"] = (Missing_gl_df["Hmax"] + Missing_gl_df["Hmin"]) / 2
Missing_gl_df["dH"] = Missing_gl_df["Hmax"] - Missing_gl_df["Hmin"]

# area in km2
Missing_gl_df["AREA_km2"] = Missing_gl_df["AREA"] * 0.000001

Missing_gl_df["AREA_LIA"] = Missing_gl_df["AREA_km2"] / 0.317

# Calculate the length
Missing_gl_df["L0_m"] = Missing_gl_df["LENGTH"]

# Calculate average slope
Missing_gl_df["a"] = Missing_gl_df.apply(
    lambda row: math.degrees(math.atan(row["dH"] / row["L0_m"])), axis=1
)

# Calculate the basal shear stress
# threshold
Missing_gl_df["tf"] = (
    0.005
    + 1.598 * (Missing_gl_df["dH"] / 1000)
    - 0.435 * ((Missing_gl_df["dH"] / 1000) ** 2)
)
Missing_gl_df.loc[Missing_gl_df["dH"] > 1600, "tf"] = 1.5
Missing_gl_df["tf_Pa"] = Missing_gl_df["tf"] * 10 ** 5

# Calculate sinus alpha
def sine(x):
    return math.sin(math.radians(x))


Missing_gl_df["sin_a"] = Missing_gl_df["a"].apply(sine)

# Calculate average ice depth along flowline (with or without f??)
Missing_gl_df["hf"] = Missing_gl_df["tf_Pa"] / (f * p * g * Missing_gl_df["sin_a"])

# Calculate average thickness
Missing_gl_df["hF"] = (math.pi / 4) * Missing_gl_df["hf"]

# Calculate maximum thickness
Missing_gl_df["h_max"] = 2.5 * Missing_gl_df["hf"]

# Calculate Volume
Missing_gl_df["V"] = Missing_gl_df["AREA_km2"] * 1000000 * Missing_gl_df["hF"]
Missing_gl_df["V_km3"] = Missing_gl_df["V"] / 1e9
Missing_gl_df["V_LIA"] = Missing_gl_df["AREA_LIA"] * 1000000 * Missing_gl_df["hF"]
Missing_gl_df["V_LIA_km3"] = Missing_gl_df["V_LIA"] / 1e9

SUM_missing_gl_LIA = Missing_gl_df["V_LIA_km3"].sum()
SUM_missing_gl = Missing_gl_df["V_km3"].sum()
SUM_missing_area = Missing_gl_df["AREA_LIA"].sum()

print(SUM_missing_gl_LIA)
print(SUM_missing_gl)
print(SUM_missing_area)



# %% Part II: Calculation of glacier volume and area changes
# %% Calcualte topographic changes, mass balance changes and annual thickness change
Combined_df = pd.merge(
    pd.merge(LIA_df, RGI_df, on="LIA_ID", suffixes=["_LIA", "_RGI"]),
    cci_df.rename(lambda x: x + "_cci" if x != "LIA_ID" else x, axis=1),
    on="LIA_ID",
)

# mass balance disturbance
Combined_df["delta_b"] = (
    Combined_df["bt_LIA"]
    * (Combined_df["L0_cci"] - Combined_df["L0_LIA"])
    / Combined_df["L0_LIA"]
)

# annual ice thickness change and conversion from m w.e to meters 
Combined_df["<b>"] = ((Combined_df["delta_b"] / 2 * Combined_df["t_resp_LIA"]) / 165)/0.9
#Combined_df["<b>"] = (Combined_df["delta_b"] / 2 * 1)/0.9

Combined_df["delta_ELA"] = Combined_df["Hmean_cci"] - Combined_df["Hmean_LIA"]
Combined_df["delta_T"] = Combined_df["delta_ELA"] / 140
mean_d_ELA = Combined_df["delta_ELA"].mean()
mean_d_T = Combined_df["delta_T"].mean()

print("ELA change"), print(mean_d_ELA)
print("mean temterature change"), print(mean_d_T)


#Mean elevation change (from total volumes) per region and both methods
region_summary_combined = pd.merge(region_summary_LIA, region_summary_cci, on="region_id", suffixes=["_LIA", "_cci"])
region_summary_combined["dV"]= region_summary_combined["V_km3_cci"]-region_summary_combined["V_km3_LIA"]
region_summary_combined["A2"]= (region_summary_combined["AREA_cci"]+region_summary_combined["AREA_LIA"])/2
region_summary_combined["dhdt"]= region_summary_combined["dV"]/region_summary_combined["A2"]/(2015-1850)*1000
dhdt_alps = (region_summary_combined["dV"].sum())/(region_summary_combined["A2"].sum())/(2015-1850)*1000
print(dhdt_alps)

## Change in mean slope
Combined_df["delta_a"] = Combined_df["a_LIA"]-Combined_df["a_cci"]

mean_slope_LIA = Combined_df["a_LIA"].mean()
mean_slpope_2015 = Combined_df["a_cci"].mean()
mean_d_slope= Combined_df["delta_a"].mean()
print("Mean slope LIA"), print(mean_slope_LIA)
print("Mean slope 2015"), print(mean_slpope_2015)
print("Slope change"), print(mean_d_slope)


# Convert to integers
Combined_df["a_LIA"] = Combined_df["a_LIA"].astype(int)
Combined_df["a_cci"] = Combined_df["a_cci"].astype(int)

mean_slope_LIA = Combined_df["a_LIA"]
mean_slope_2015 = Combined_df["a_cci"]

# Fit a linear regression line
slope, intercept, r_value, p_value, std_err = linregress(mean_slope_LIA, mean_slope_2015)
trendline = slope * mean_slope_2015 + intercept

# Plot the scatter plot and the linear trendline
plt.figure()
plt.scatter(mean_slope_LIA, mean_slope_2015, 
            marker="o", 
            edgecolors="black", 
            facecolors="none",
            alpha= 0.5, 
            label = "glaciers; " f"Linear Trendline: y = {slope:.2f} * x + {intercept:.2f}, R² = {r_value:.2f} ")
plt.plot(mean_slope_2015, 
         trendline, 
         color='red')

plt.xlabel("Mean slope LIA (°)")
plt.ylabel("Mean slope 2015 (°)")
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.legend()
plt.show()

#Length changes
Combined_df["delta_L"] = Combined_df["L0_cci"]-Combined_df["L0_LIA"]


# %% Volume changes
Volume_change = cci_volume - LIA_volume
print('Volume change'), print(Volume_change)

Volume_change_regions = cci_sum_vol_region- LIA_sum_vol_region
print("\Volume change per region:")
print(Volume_change_regions)

# %% Area changes
LIA_df["AREA"].replace([np.inf, -np.inf, np.nan], 0, inplace=True)
LIA_area = LIA_df["AREA"].sum()
RGI_area = RGI_df.loc[RGI_df["LIA_ID"] < 5000, "AREA"].sum()
RGI_area_all = RGI_df["AREA"].sum()
cci_area = cci_df.loc[cci_df["LIA_ID"] < 5000, "AREA"].sum()
cci_area_all = cci_df["AREA"].sum()

print("LIA ice area"), print(LIA_area)
print("RGI ice area"), print(RGI_area)
print("RGI ice area_all"), print(RGI_area_all)
print("cci ice area"), print(cci_area)
print("cci ice area_all"), print(cci_area_all)
Area_change = cci_area - LIA_area
print("Area change"), print(Area_change)

# Calculate relative changes for each region between LIA and RGI
relative_change_LIA_RGI = (100 - (100 / LIA_area * RGI_area)) * -1
relative_change_LIA_cci = (100 - (100 / LIA_area * cci_area)) * -1
relative_change_RGI_cci = relative_change_LIA_cci - relative_change_LIA_RGI
relative_change_rate_LIA_cci = relative_change_LIA_cci / (2015 - 1850)
relative_change_rate_LIA_RGI = relative_change_LIA_RGI / (2003 - 1850)
relative_change_rate_RGI_cci = (100 - (100 / RGI_area * cci_area)) * -1 / (2015 - 2003)

print("Relative area change LIA-RGI"), print(relative_change_LIA_RGI)
print("Relative area change LIA-cci"), print(relative_change_LIA_cci)
print("Relative area change RGI-cci"), print(relative_change_RGI_cci)
print("Relative area change rate LIA-cci"), print(relative_change_rate_LIA_cci)
print("Relative area change rate LIA-RGI"), print(relative_change_rate_LIA_RGI)
print("Relative area change rate RGI-cci"), print(relative_change_rate_RGI_cci)


# Area change per region
# Step 1: Calculate the sum of "AREA" per region for each dataframe
LIA_sum = LIA_df.groupby("region_id")["AREA"].sum().reset_index()
RGI_sum = RGI_df.groupby("region_id")["AREA"].sum().reset_index()
cci_sum = cci_df.groupby("region_id")["AREA"].sum().reset_index()
sxs_sum = SxS_df.groupby("region_id")["AREA"].sum().reset_index()

# Step 2: Merge the dataframes and calculate relative changes compared to LIA_df
area_region = LIA_sum.merge(
    RGI_sum, on="region_id", suffixes=("_LIA", "_RGI"), how="outer"
)
area_region = area_region.merge(cci_sum, on="region_id", how="outer")

# Calculate relative changes compared to LIA_df
area_region["Relative area change LIA-RGI"] = (
    100 - (100 / area_region["AREA_LIA"] * area_region["AREA_RGI"])
) * -1
area_region["Relative area change LIA-CCI"] = (
    100 - (100 / area_region["AREA_LIA"] * area_region["AREA"])
) * -1
area_region["Relative area change RGI-CCI"] = (
    area_region["Relative area change LIA-CCI"]
    - area_region["Relative area change LIA-RGI"]
)

area_region["region"] = area_region["region_id"].replace(
    {
        1: "west",
        2: "west",
        3: "west",
        4: "west",
        5: "west",
        6: "west",
        7: "west",
        8: "west",
        9: "east",
        10: "east",
        11: "east",
        12: "east",
        13: "east",
        14: "east",
    }
)

# Now you can calculate the sums for 'west' and 'east'
sum_area_west = area_region[area_region["region"] == "west"].sum()
sum_area_east = area_region[area_region["region"] == "east"].sum()

# Extract the sum values for each AREA column
sum_lia_west = sum_area_west["AREA_LIA"]
sum_rgi_west = sum_area_west["AREA_RGI"]
sum_area_west = sum_area_west["AREA"]

sum_lia_east = sum_area_east["AREA_LIA"]
sum_rgi_east = sum_area_east["AREA_RGI"]
sum_area_east = sum_area_east["AREA"]
area_change_west = (100 - (100 / sum_lia_west * sum_area_west)) * -1
area_change_east = (100 - (100 / sum_lia_east * sum_area_east)) * -1

area_change_west_RGI = (100 - (100 / sum_lia_west * sum_rgi_west)) * -1
area_change_east_RGI = (100 - (100 / sum_lia_east * sum_rgi_east)) * -1

area_change_west_cci = area_change_west - area_change_west_RGI
area_change_east_cci = area_change_east - area_change_east_RGI

area_change_rate_west_LIA_RGI = area_change_west_RGI / (2003 - 1850)
area_change_rate_east_LIA_RGI = area_change_east_RGI / (2003 - 1850)
area_change_rate_west_RGI_cci = ((100 - (100 / sum_rgi_west * sum_area_west)) * -1) / (
    2015 - 2003
)
area_change_rate_east_RGI_cci = ((100 - (100 / sum_rgi_east * sum_area_east)) * -1) / (
    2015 - 2003
)


# Print the sums for 'west' and 'east'
print("Sum of AREA_LIA for 'west':", sum_lia_west)
print("Sum of AREA_RGI for 'west':", sum_rgi_west)
print("Sum of AREA for 'west':", sum_area_west)
print("Sum of AREA_LIA for 'east':", sum_lia_east)
print("Sum of AREA_RGI for 'east':", sum_rgi_east)
print("Sum of AREA for 'east':", sum_area_east)
print("area_change_west':", area_change_west)
print("area_change_east':", area_change_east)
print("area_change_west_LIA_RGI':", area_change_west_RGI)
print("area_change_east_LIA_RGI':", area_change_east_RGI)
print("area_change_west_RGI_cci':", area_change_west_cci)
print("area_change_east_RGI_cci':", area_change_east_cci)
print("area_change_rate_west_LIA_RGI':", area_change_rate_west_LIA_RGI)
print("area_change_rate_east_LIA_RGI':", area_change_rate_east_LIA_RGI)
print("area_change_rate_west_RGI_cci':", area_change_rate_west_RGI_cci)
print("area_change_rate_east_RGI_cci':", area_change_rate_east_RGI_cci)

area_region["LIA_100"] = 0
area_region["Date_LIA"] = 1850
area_region["Date_RGI"] = 2003
area_region["Date_cci"] = 2015

area_region["rate_LIA_RGI"] = area_region["Relative area change LIA-RGI"] / (
    area_region["Date_RGI"] - area_region["Date_LIA"]
)
area_region["rate_RGI_cci"] = (
    (100 - (100 / area_region["AREA_RGI"] * area_region["AREA"])) * -1
) / (area_region["Date_cci"] - area_region["Date_RGI"])

# Area change size classes cci
merged_df_AREA_cci = pd.merge(LIA_df, cci_df, on="LIA_ID", suffixes=("_LIA", "_cci"))
merged_df_AREA_cci["RELATIVE_AREA_CHANGE"] = (
    100 - (100 / merged_df_AREA_cci["AREA_LIA"] * merged_df_AREA_cci["AREA_cci"])
) * -1
# Create bins for different size classes
bins = [0, 1, 5, 10, 50, float("inf")]
labels = ["<1", "1-5", "5-10", "10-50", ">50"]
# Cut the data into different size classes based on LIA_df['SIZE_CLASS']
merged_df_AREA_cci["LIA_SIZE_CLASS"] = pd.cut(
    merged_df_AREA_cci["AREA_LIA"], bins=bins, labels=labels
)
# Calculate relative area change for each size class
sum_by_size_class = merged_df_AREA_cci.groupby("LIA_SIZE_CLASS")[
    ["AREA_LIA", "AREA_cci"]
].sum()
# Calculate relative change for each size class
sum_by_size_class["RELATIVE_AREA_CHANGE"] = (
    100 - (100 / sum_by_size_class["AREA_LIA"] * sum_by_size_class["AREA_cci"])
) * -1
print(sum_by_size_class)

# Area change size classes RGI
merged_df_AREA_RGI = pd.merge(LIA_df, RGI_df, on="LIA_ID", suffixes=("_LIA", "_RGI"))
merged_df_AREA_RGI["RELATIVE_AREA_CHANGE"] = (
    100 - (100 / merged_df_AREA_RGI["AREA_LIA"] * merged_df_AREA_RGI["AREA_RGI"])
) * -1
# Cut the data into different size classes based on LIA_df['SIZE_CLASS']
merged_df_AREA_RGI["LIA_SIZE_CLASS"] = pd.cut(
    merged_df_AREA_RGI["AREA_LIA"], bins=bins, labels=labels
)
# Calculate relative area change for each size class
sum_by_size_class2 = merged_df_AREA_RGI.groupby("LIA_SIZE_CLASS")[
    ["AREA_LIA", "AREA_RGI"]
].sum()
# Calculate relative change for each size class
sum_by_size_class2["RELATIVE_AREA_CHANGE"] = (
    100 - (100 / sum_by_size_class2["AREA_LIA"] * sum_by_size_class2["AREA_RGI"])
) * -1
print(sum_by_size_class2)

# ------------------------------------------------------------------------------
# calculate area changes and area change rates
Combined_df["Area_change_rate_LIA_RGI"] = (
    (100 - (100 / Combined_df["AREA_LIA"] * Combined_df["AREA_RGI"])) * -1
) / (2003 - 1850)
Combined_df["Area_change_rate_RGI_cci"] = (
    (100 - (100 / Combined_df["AREA_RGI"] * Combined_df["AREA_cci"])) * -1
) / (2015 - 2003)
Combined_df["Area_change_rate_LIA_cci"] = (
    (100 - (100 / Combined_df["AREA_LIA"] * Combined_df["AREA_cci"])) * -1
) / (2015 - 1850)


Combined_df["Area_change_LIA_RGI"] = (
    100 - (100 / Combined_df["AREA_LIA"] * Combined_df["AREA_RGI"])
) * -1
Combined_df["Area_change_LIA_cci"] = (
    100 - (100 / Combined_df["AREA_LIA"] * Combined_df["AREA_cci"])
) * -1
Combined_df["Area_change_RGI_cci"] = (
    100 - (100 / Combined_df["AREA_RGI"] * Combined_df["AREA_cci"])
) * -1

# ------------------------------------------------------------------------------
# glacier specific area change and std

print("glacier specific area change RGI-cci:", Combined_df["Area_change_RGI_cci"].mean(),"+-",Combined_df["Area_change_RGI_cci"].std())

# %% Calculate sea level rise equivalent
LIA_df["Gt"] = LIA_df["V_km3"] * 0.9
RGI_df["Gt"] = RGI_df["V_km3"] * 0.9

# Calculate SLE
LIA_df["Gt"].replace([np.inf, -np.inf, np.nan], 0, inplace=True)
RGI_df["Gt"].replace([np.inf, -np.inf, np.nan], 0, inplace=True)
SLE_past = (LIA_df["Gt"].sum() * (1 / 361.8)) - (RGI_df["Gt"].sum()) * (1 / 361.8)
SLE_total = LIA_df["Gt"].sum() * (1 / 361.8)
SLE_potential = RGI_df["Gt"].sum() * (1 / 361.8)
print("SLE_past"), print(SLE_past)
print("SLE_total"), print(SLE_total)
print("SLE_potential"), print(SLE_potential)





# Part III: plots
# %% Area change plots
## plots change per size class
sum_by_size_class_merged = pd.merge(
    sum_by_size_class,
    sum_by_size_class2,
    on="LIA_SIZE_CLASS",
    suffixes=("_cci", "_RGI"),
)

# Assuming the index is used as the x-axis
x_axis_column = sum_by_size_class_merged.index
y_axis_columns = ["RELATIVE_AREA_CHANGE_cci", "RELATIVE_AREA_CHANGE_RGI"]

# Extract relevant data from the DataFrame
y_values = sum_by_size_class_merged[y_axis_columns]

# Plotting using seaborn
plt.figure(figsize=(10, 6))
sns.barplot(
    x=x_axis_column,
    y=y_values[y_axis_columns[0]],
    data=sum_by_size_class_merged,
    color="black",
    label="2003-2015",
)
sns.barplot(
    x=x_axis_column,
    y=y_values[y_axis_columns[1]],
    data=sum_by_size_class_merged,
    color="grey",
    label="LIA-2003",
)

for x in range(len(x_axis_column) - 1):
    plt.axvline(x + 0.5, color="black", linestyle="--", linewidth=0.8)

plt.grid(axis="y", linestyle="--", linewidth=0.5)

# Adding labels and title
plt.xlabel("Size class (km$^{2}$)")
plt.ylabel("Relative area change (%)")

# Adding legend
plt.legend(title="Legend")

# Show the plot
plt.show()

# ------------------------------------------------------------------------------
# plot area changes per region
colors  = [
    "#FF0000",
    "#FF6A00",
    "#FFD400",
    "#C0FF00",
    "#56FF00",
    "#00FF14",
    "#00FF7E",
    "#00FFE7",
    "#00ADFF",
    "#0043FF",
    "#2700FF",
    "#9100FF",
    "#FB00FF",
    "#FF0099"
]

shape_mapping = {
    "1": "s",       # Filled triangle for "1"
    "2": "o",       # Filled plus sign for "2"
    "3": "^",       # Filled cross for "3"
    "4": "D",       # Filled diamond for "4"
    "5": "s",       # Filled square for "5"
    "6": "o",       # Filled octagon for "6"
    "7": "^",       # Filled star for "7"
    "8": "D",       # Filled circle for "8"
    "9": "s",       # Filled hexagon for "9"
    "10": "o",      # Filled circle for "10"
    "11": "^",      # Empty triangle for "11"
    "12": "D",      # Empty plus sign for "12"
    "13": "s",      # Empty cross for "13"
    "14": "o",      # Empty star for "14"
}

shape_fill_styles = {
    "1": "none",    # Empty for "1"
    "2": "full",    # Filled for "2"
    "3": "none",    # Filled for "3"
    "4": "full",    # Filled for "4"
    "5": "none",    # Filled for "5"
    "6": "full",    # Filled for "6"
    "7": "none",    # Filled for "7"
    "8": "none",    # Filled for "8"
    "9": "full",    # Filled for "9"
    "10": "none",   # Filled for "10"
    "11": "full",   # Empty for "11"
    "12": "none",   # Empty for "12"
    "13": "full",   # Empty for "13"
    "14": "none",   # Empty for "14"
}


# Sort the dataframe by region_id for better visualization of connecting lines
area_region.sort_values(by="region_id", inplace=True)

# Create a list of unique region IDs for iteration
region_ids = area_region["region_id"].unique()

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 6))
x_ticks = np.arange(1850, 2026, 25)

# Iterate through each region and plot the points and connecting lines
for i, region_id in enumerate(region_ids):
    region_data = area_region[area_region["region_id"] == region_id]
    x = [
        region_data["Date_LIA"].values[0],
        region_data["Date_RGI"].values[0],
        region_data["Date_cci"].values[0],
    ]
    y = [
        region_data["LIA_100"].values[0],
        region_data["Relative area change LIA-RGI"].values[0],
        region_data["Relative area change LIA-CCI"].values[0],
    ]

    # Determine marker style based on region_id
    marker = shape_mapping.get(str(region_id), 'o')

    # Retrieve fill style based on region_id
    fill_style = shape_fill_styles.get(str(region_id), 'None')

    # Plot the points with the specified marker and fill style
    # Plot the points with the specified marker and fill style
    for j, (x_val, y_val) in enumerate(zip(x, y)):
        if fill_style == 'none':
            ax.plot(x_val, y_val, marker=marker, color=colors[i % len(colors)], markersize=7, markeredgewidth=1.5, linestyle='None', markerfacecolor='none', markeredgecolor=colors[i % len(colors)])
        else:
            ax.plot(x_val, y_val, marker=marker, color=colors[i % len(colors)], markersize=7, linestyle='None')
    
    # Plot the connecting lines
    ax.plot(x, y, linestyle="-", marker=marker, color=colors[i % len(colors)], label=f"Region {region_id}")


# Set labels, ticks, legend, and grid
ax.set_xlabel("Year")
ax.set_ylabel("Relative area change (%)")
ax.set_xticks(x_ticks)
ax.set_xticklabels([str(year) for year in x_ticks])
ax.legend()
plt.grid(True)

# Show the plot
plt.show()

# ------------------------------------------------------------------------------
# #plot area change rates per region

colors  = [
    "#FF0000",
    "#FF6A00",
    "#FFD400",
    "#C0FF00",
    "#56FF00",
    "#00FF14",
    "#00FF7E",
    "#00FFE7",
    "#00ADFF",
    "#0043FF",
    "#2700FF",
    "#9100FF",
    "#FB00FF",
    "#FF0099"
]

shape_mapping = {
    "1": "s",       # Filled triangle for "1"
    "2": "o",       # Filled plus sign for "2"
    "3": "^",       # Filled cross for "3"
    "4": "D",       # Filled diamond for "4"
    "5": "s",       # Filled square for "5"
    "6": "o",       # Filled octagon for "6"
    "7": "^",       # Filled star for "7"
    "8": "D",       # Filled circle for "8"
    "9": "s",       # Filled hexagon for "9"
    "10": "o",      # Filled circle for "10"
    "11": "^",      # Empty triangle for "11"
    "12": "D",      # Empty plus sign for "12"
    "13": "s",      # Empty cross for "13"
    "14": "o",      # Empty star for "14"
}
area_region.sort_values(by="region_id", inplace=True)
region_ids = area_region["region_id"].unique()

fig, ax = plt.subplots(figsize=(10, 6))
x_ticks = np.arange(1850, 2026, 25)

custom_lines = []

for i, region_id in enumerate(region_ids):
    region_data = area_region[area_region["region_id"] == region_id]

    x1 = region_data["Date_LIA"].values[0]
    x2 = region_data["Date_RGI"].values[0]
    x3 = region_data["Date_cci"].values[0]

    y1 = region_data["rate_LIA_RGI"].values[0]
    y2 = region_data["rate_RGI_cci"].values[0]

    marker = shape_mapping.get(str(region_id), 'o')


    # Plot the points with the defined color palette and markers
    ax.step(
        [x1, x2, None, x2, x3],
        [y1, y1, None, y2, y2],
        where="post",
        color=colors[i % len(colors)],
        label=f"Region {region_id}",
    )
    ax.scatter(
        [x1, x2, x2, x3],
        [y1, y1, y2, y2],
        marker=marker,
        s=50,
        color=colors[i % len(colors)],
        zorder=5,
        label=f"Region {region_id}",
    )

    # Create custom legend handler
    line = Line2D(
        [0],
        [0],
        color=colors[i % len(colors)],
        label=f"Region {region_id}",
        linestyle="-",
        marker=marker,
    )
    custom_lines.append(line)

# Add custom legend entries for lines and markers
ax.legend(handles=custom_lines, loc="lower left")

ax.set_xlabel("Year")
ax.set_ylabel("Relative area change rate (m a$^{-1}$)")
ax.set_xticks(x_ticks)
ax.set_xticklabels([str(year) for year in x_ticks])
plt.grid(True)
plt.show()

# ------------------------------------------------------------------------------
# scatter plot of relative area change rate LIA-2015 against glacier size
plt.figure(figsize=(8, 6))
plt.scatter(
    Combined_df["AREA_LIA"],
    Combined_df["Area_change_rate_LIA_cci"],
    marker="o",
    facecolors="none",
    edgecolors="black",
)
plt.xscale("log")
# Set axis labels and title
plt.xlabel("Area LIA (km$^{2}$)")
plt.ylabel("Relative area change rate (% a$^{-1}$)")
plt.axhline(0, color="black", linestyle="-", linewidth=1)
plt.ylim(-0.7, 0.2)
# Show the plot
plt.grid(True)
plt.show()

# ------------------------------------------------------------------------------
# scatter plot of relative area change LIA-2015 against glacier size
plt.figure(figsize=(8, 6))
plt.scatter(
    Combined_df["AREA_LIA"],
    Combined_df["Area_change_LIA_cci"],
    marker="o",
    facecolors="none",
    edgecolors="black",
)
plt.xscale("log")
# Set axis labels and title
plt.xlabel("Area LIA (km$^{2}$)")
plt.ylabel("Relative area change (%)")
plt.axhline(0, color="black", linestyle="-", linewidth=1)
plt.ylim(-100, 20)
# Show the plot
plt.grid(True)
plt.show()

# ------------------------------------------------------------------------------
# plot area change rates against glacier size for both time periods
plt.figure(figsize=(8, 6))
plt.scatter(
    Combined_df["AREA_RGI"],
    Combined_df["Area_change_rate_RGI_cci"],
    marker="o",
    color="black",
)
plt.scatter(
    Combined_df["AREA_LIA"],
    Combined_df["Area_change_rate_LIA_RGI"],
    marker="o",
    facecolors="none",
    edgecolors="red",
)
plt.xscale("log")
# Set axis labels and title
plt.xlabel("Area (km$^{2}$)")
plt.ylabel("Relative area change rate (% a$^{-1}$)")
plt.ylim(-7, 2)
plt.legend(["2003-2015", "LIA-2003"])

# Show the plot
plt.grid(True)
plt.show()


# %% Volume and elevation change plots ------------------------------------------------------------------------------
combined_df_methods = pd.merge(Parameter_df, GIS_df, on="LIA_ID", how="inner")
# filtered_columns = [col for col in combined_df_methods.columns if col.endswith('_GIS')]
# combined_df_methods = pd.merge(Combined_df, combined_df_methods[['LIA_ID']], on='LIA_ID', how='inner')
combined_df_methods["Difference_methods"] = (
    combined_df_methods["delta_V"] - combined_df_methods["Volume_change_GIS"]
)


# Plot the difference between the methods against the slope
combined_df_methods["slope"] = pd.to_numeric(
    combined_df_methods["slope"], errors="coerce"
)
combined_df_methods["Difference_methods"] = pd.to_numeric(
    combined_df_methods["Difference_methods"], errors="coerce"
)
filtered_df = combined_df_methods.dropna(subset=["slope", "Difference_methods"])

# Perform linear regression
slope_mean = np.mean(filtered_df["slope"])
difference_mean = np.mean(filtered_df["Difference_methods"])
slope_diff = filtered_df["slope"] - slope_mean
difference_diff = filtered_df["Difference_methods"] - difference_mean
slope_diff_squared = slope_diff * slope_diff
slope_diff_difference_diff = slope_diff * difference_diff
slope_regression = slope_diff_difference_diff.sum() / slope_diff_squared.sum()
intercept_regression = difference_mean - slope_regression * slope_mean

plt.figure(figsize=(8, 6))
plt.scatter(filtered_df["slope"], filtered_df["Difference_methods"],  s=5, c='black', alpha=0.5)
plt.xlabel("Mean slope (°)", fontsize = 17)
plt.ylabel("Difference between methods (km$^3$)",fontsize = 17)
ax.tick_params(axis='both', which='major', labelsize=18)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.xlim(0, 60)
plt.ylim(-0.5, 0.5)
plt.grid(True)

# Add the linear regression line
plt.plot(
    filtered_df["slope"],
    slope_regression * filtered_df["slope"] + intercept_regression,
    color="grey",
)
plt.show()


# --------------------------------------------------------------------------------------------------------------------
# Plot the difference between the methods against the AREA
combined_df_methods["AREA_LIA"] = pd.to_numeric(
    combined_df_methods["AREA_LIA"], errors="coerce"
)
combined_df_methods["Difference_methods"] = pd.to_numeric(
    combined_df_methods["Difference_methods"], errors="coerce"
)
filtered_df = combined_df_methods.dropna(subset=["AREA_LIA", "Difference_methods"])

# Perform linear regression
area_mean = np.mean(filtered_df["AREA_LIA"])
area_diff = filtered_df["AREA_LIA"] - area_mean
area_diff_squared = area_diff * area_diff
area_diff_difference_diff = area_diff * difference_diff
area_regression = area_diff_difference_diff.sum() / area_diff_squared.sum()
intercept_regression = difference_mean - area_regression * area_mean

plt.figure(figsize=(8, 6))
plt.scatter(filtered_df["AREA_LIA"], filtered_df["Difference_methods"], s=5, c='black',alpha=0.5)
plt.xlabel("Area during LIA (km$^2$)",fontsize = 17)
plt.ylabel("Difference between methods (km$^3$)",fontsize = 17)
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.xlim(0, 20)
plt.ylim(-0.5, 0.5)
plt.grid(True)
# Add the linear regression line
plt.plot(
    filtered_df["AREA_LIA"],
    area_regression * filtered_df["AREA_LIA"] + intercept_regression,
    color="grey", linestyle="--",
)
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# Plot Volume change per glacier for both methods
x = combined_df_methods["Volume_change_GIS"]
y = combined_df_methods["delta_V"]

# x = x.drop(1895)
# y = y.drop(1895)

# Calculate the linear trendline and  R-squared value
coefficients = np.polyfit(x, y, 1)
a, b = coefficients
y_predicted = a * x + b
residuals = y - y_predicted
ss_residuals = np.sum(residuals ** 2)
ss_total = np.sum((y - np.mean(y)) ** 2)
r_squared = 1 - (ss_residuals / ss_total)

# Scatter plot with linear scales
plt.figure(figsize=(8, 6))
plt.scatter(x, y, alpha=0.5,marker="o", edgecolors="black", facecolors="none", s= 10)
plt.xlim(-1, 0)
plt.ylim(-1, 0)
plt.xlabel("Volume change GIS (km$^3$)", fontsize = 14)
plt.ylabel("Volume change parameter (km$^3$)", fontsize = 14)
x_range = np.linspace(min(x), max(x), 100)
y_trendline = a * x_range + b
plt.plot(
    x_range,
    y_trendline,
    color="red",
    label=f"Linear Trendline (y = {a:.2f} * x + {b:.2f}, R² = {r_squared:.2f} )",
)
plt.plot([-7, 0], [-7, 0], "k-", linestyle="--")
plt.legend()
plt.grid(True)
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# plot of mean elevation change rates per region and method against mean elevation
# Calculate mean values per region and method
mean_values = combined_df_methods.groupby(["reg_ID_NEW"]).mean().reset_index()

# Define your color palette
# colors = ['#FF6A00', '#FFD400', '#C0FF00', '#56FF00', '#00FF14', '#00FF7E', '#00FFE7', '#00ADFF', '#0043FF', '#0043FF', '#2700FF', '#9100FF', '#FB00FF', '#FF0099']
gist_rainbow_colors = plt.cm.gist_rainbow(np.linspace(0, 1, 14))

# Converting RGBA to hex
colors = [
    "#%02x%02x%02x" % tuple(int(c * 255) for c in rgba[:-1])
    for rgba in gist_rainbow_colors
]

# Create a scatter plot for mean values
for i, region_id in enumerate(mean_values["reg_ID_NEW"].unique()):
    region_data = mean_values[mean_values["reg_ID_NEW"] == region_id]
    color = colors[i % len(colors)]  # Use the color palette
    plt.scatter(
        region_data["<b>"],
        region_data["MEAN_LIA"],
        marker="x" if region_id in range(1, 9) else "+",
        label=f"Parameter approach - Region {region_id}",
        alpha=1,
        color=color,
        s=60,
    )

# Create a separate scatter plot for the 'GIS approach'
for i, region_id in enumerate(mean_values["reg_ID_NEW"].unique()):
    region_data = mean_values[mean_values["reg_ID_NEW"] == region_id]
    color = colors[i % len(colors)]  # Use the color palette
    plt.scatter(
        region_data["dhdt_GIS"],
        region_data["MEAN_LIA"],
        marker="o" if region_id in range(1, 9) else "s",
        alpha=1,
        color=color,
        s=60,
    )

legend_elements = [
    Line2D(
        [0],
        [0],
        marker="x",
        color="black",
        markersize=10,
        label="Parameter approach - Western Alps",
    ),
    Line2D(
        [0],
        [0],
        marker="+",
        color="black",
        markersize=10,
        label="Parameter approach - Eastern Alps",
    ),
    Line2D(
        [0],
        [0],
        marker="o",
        color="black",
        markersize=10,
        label="GIS approach - Western Alps",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="black",
        markersize=10,
        label="GIS approach - Eastern Alps",
    ),
]
plt.gca().add_artist(
    plt.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=True,
        facecolor="w",
        edgecolor="w",
        framealpha=0.4,
    )
)


# Custom legend for regions 1-14 with specific colors represented by lines
region_legend_elements = [
    Line2D([0], [0], color=colors[i % len(colors)], lw=2, label=f" {i+1}")
    for i in range(14)
]
plt.gca().add_artist(
    plt.legend(
        handles=region_legend_elements,
        title="Regions",
        loc="upper left",
        frameon=True,
        facecolor="w",
        edgecolor="w",
        framealpha=0.4,
    )
)

plt.xlabel("Mean annual thickness change (m a$^{-1}$)")
plt.ylabel("Mean elevation (m)")
plt.xlim(-0.3, 0)
plt.ylim(2000, 3200)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# plot Mass balance at the tongue agains min elevation
plt.figure()
plt.scatter(
    LIA_df["bt"], LIA_df["MIN"], marker="o", edgecolors="black", alpha = 0.5, facecolors="none",s= 10
)
plt.xlabel("Mass balance at glacier terminus (m a$^{-1}$)", fontsize = 14)
plt.ylabel("Minimum elevation (m)", fontsize = 14)
x = LIA_df["bt"]
y = LIA_df["MIN"]
slope, intercept, r_value, p_value, std_err = linregress(x, y)
plt.plot(
    x,
    slope * x + intercept,
    color="red",
    linestyle="--",
    label=f"Linear Regression (R-squared={r_value**2:.2f})",
)
plt.legend()
plt.show()

# plot histogram of mass balance at glacier tongue
plt.figure()
plt.hist(LIA_df['bt'], bins=20, edgecolor='black', alpha=0.7)
plt.xlabel('Mass balance at glacier tongue (m a$^{-1}$)')
plt.ylabel('Frequency')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()
# --------------------------------------------------------------------------------------------------------------------
# Mean glacier specific annual thickness change per region and method
mean_dhdt_param = Parameter_df["<b>"].mean()
print(mean_dhdt_param)
mean_dhdt_GIS = GIS_df["dhdt_GIS"].mean()
print(mean_dhdt_GIS)

mean_values = Parameter_df.groupby("region_id_cci")["<b>"].mean()
std_values = Parameter_df.groupby("region_id_cci")["<b>"].std()
mean_values_gis = GIS_df.groupby("reg_ID_NEW")["dhdt_GIS"].mean()
std_values_gis = GIS_df.groupby("reg_ID_NEW")["dhdt_GIS"].std()

fig, ax = plt.subplots()
width = 0.4  # Adjust the width for bars
x = range(len(mean_values))  # Create x values for the bars
    
# Plot the Parameter data
ax.bar(
    x,
    mean_values,
    yerr=std_values,
    capsize=4,
    width=width,
    color="blue",
    alpha=0.7,
    edgecolor="black",
    label="Parameter",
)

# Plot the GIS data
ax.bar(
    [i + width for i in x],
    mean_values_gis,
    yerr=std_values_gis,
    capsize=4,
    width=width,
    color="yellow",
    alpha=0.7,
    edgecolor="black",
    label="GIS",
)

# Draw dashed lines between regions
for i in range(0, len(mean_values)):
    ax.axvline(x=i + 0.3 + width, color="black", linestyle="--", linewidth=0.5)

ax.set_xlabel("Subregions (Region ID)")
ax.set_ylabel("Mean annual thickness change (m a$^{-1}$)")
ax.legend()
plt.xticks(
    [i + width / 2 for i in x], mean_values.index
)  # Set x-ticks at the center of each pair of bars
plt.show()

#--------------------------------------------------------------------------------------------------------------------
# Mean annual thickness change per region and method
mean_values = region_summary_combined["dhdt"]
print(mean_values)
mean_values_GIS = GIS_region["dh_GIS"]
selected_values_GIS = mean_values_GIS.iloc[:14]
print(selected_values_GIS)


mean_values_gis = GIS_df.groupby("reg_ID_NEW")["dhdt_GIS"].mean()


fig, ax = plt.subplots()
width = 0.4  # Adjust the width for bars
x = range(len(mean_values))  # Create x values for the bars
    
# Plot the Parameter data
ax.bar(
    x,
    mean_values,
    capsize=4,
    width=width,
    color="blue",
    alpha=0.7,
    edgecolor="black",
    label="Parameter",
)

# Plot the GIS data
ax.bar(
    [i + width for i in x],
    selected_values_GIS,
    capsize=4,
    width=width,
    color="yellow",
    alpha=0.7,
    edgecolor="black",
    label="GIS",
)

# Draw dashed lines between regions
for i in range(0, len(mean_values)):
    ax.axvline(x=i + 0.3 + width, color="black", linestyle="--", linewidth=0.5)

# Add dashed horizontal line of Alpine means
dh_GIS_total = GIS_region.loc[GIS_region["region_id"] == "Total", "dh_GIS"].values[0]
ax.axhline(y=dh_GIS_total, color="yellow", linestyle="--", linewidth=2)

dhdt_alps_value = dhdt_alps  # Replace 'dhdt_alps' with the actual value
ax.axhline(y=dhdt_alps_value, color="blue", linestyle="--", linewidth=2)


ax.set_xlabel("Region ID", fontsize=14)
ax.set_ylabel("Mean elevation change rate (m a$^{-1}$)", fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.legend(fontsize=12)
plt.xticks([i + width / 2 for i in x], range(1, 15),fontsize=12)
plt.show()

#optional export of the figure
#fig.set_size_inches(8.27, 6)
#plt.savefig('filepath/filename.jpeg', dpi=300)
# --------------------------------------------------------------------------------------------------------------------
# plot mean mass balance disturbance
mean_values = Combined_df.groupby("region_id_cci")["delta_b"].mean()
std_values = Combined_df.groupby("region_id_cci")["delta_b"].std()
fig, ax = plt.subplots()
mean_values.plot(
    kind="bar",
    ax=ax,
    yerr=std_values,
    capsize=4,
    color="blue",
    alpha=0.7,
    edgecolor="black",
)
ax.set_xlabel("Subregions", fontsize = 14)
ax.set_ylabel("Mean mass balance disturbance (m)", fontsize = 14)
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# bar plot area per region
lia_sum_area = LIA_df.groupby("region_id")["AREA"].sum()
RGI_sum_area = RGI_df.groupby("region_id")["AREA"].sum()
cci_sum_area = cci_df.groupby("region_id")["AREA"].sum()
combined_sum_df = pd.DataFrame(
    {
        "LIA_SUM_AREA": lia_sum_area,
        "RGI_SUM_AREA": RGI_sum_area,
        "CCI_SUM_AREA": cci_sum_area,
    }
)
combined_sum_df.plot(kind="bar", figsize=(10, 6))
plt.xlabel("Subregion")
plt.ylabel("Sum of area (km2)")
plt.title("Sum of area per Subregion ")
plt.legend(["LIA", "2003", "2015"])
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# plot reaction time agains average slope
plt.figure()
plt.scatter(
    cci_df["t_resp"], cci_df["a"], marker="o", edgecolors="red", facecolors="none", s=20, label="2015"
)
plt.scatter(
    LIA_df["t_resp"], LIA_df["a"], marker="o", edgecolors="black",alpha = 0.3,facecolors="none",s=20, label="LIA"
)
plt.xlabel("Response time (a)", fontsize=14)
plt.ylabel("Average slope (°)", fontsize=14)
plt.ylim(0, 70)
plt.xlim(0, 120)
plt.legend()
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# plot histogram of responce time
min_value = min(LIA_df["t_resp"].min(), cci_df["t_resp"].min())
max_value = max(LIA_df["t_resp"].max(), cci_df["t_resp"].max())
num_bins = 40

plt.figure()
plt.hist(LIA_df["t_resp"], bins=num_bins, range=(min_value, max_value), edgecolor="black", alpha=0.7, color="blue", label="LIA")
plt.hist(cci_df["t_resp"], bins=num_bins, range=(min_value, max_value), edgecolor="black", alpha=0.7, color="yellow", label="2015")
plt.xlabel("Response Time (a)", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.xlim(0, 120)
plt.legend()
plt.show()

# --------------------------------------------------------------------------------------------------------------------
# plot scatterplot of median elevation (GIS) vs mean from parameter
# Remove rows with zero values
combined_df_methods["MEDIAN_cci"] = combined_df_methods["MEDIAN_cci"].replace([np.inf, -np.inf, np.nan], 0).astype(int)
combined_df_methods["Hmean_cci"] = combined_df_methods["Hmean_cci"].replace([np.inf, -np.inf, np.nan], 0).astype(int)
combined_df_methods = combined_df_methods[(combined_df_methods["MEDIAN_cci"] != 0) & (combined_df_methods["Hmean_cci"] != 0)]

# Convert to integers
combined_df_methods["MEDIAN_cci"] = combined_df_methods["MEDIAN_cci"].astype(int)
combined_df_methods["Hmean_cci"] = combined_df_methods["Hmean_cci"].astype(int)

median_elevation = combined_df_methods["MEDIAN_cci"]
h_mean = combined_df_methods["Hmean_cci"]

# Fit a linear regression line
slope, intercept, r_value, p_value, std_err = linregress(median_elevation, h_mean)
trendline = slope * median_elevation + intercept

# Plot the scatter plot and the linear trendline
plt.figure()
plt.scatter(median_elevation, h_mean, marker="o", edgecolors="black", facecolors="none",alpha= 0.5, s = 10, label = "glaciers; " f" R² = {r_value:.2f} ")
plt.plot(median_elevation, 
         trendline, 
         color='red')

plt.xlabel("Median elevation (m)", fontsize = 14)
plt.ylabel("Mid-point elevation (m)", fontsize = 14)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.legend()
plt.show()


# --------------------------------------------------------------------------------------------------------------------
#scatter plot of flow velocities from Parameter scheme and data from 2015/2016 from Rabatel et al. 2023
# Assuming "LIA_ID" is the common column for merging
merged_df_flow_velocity = pd.merge(flow_df, cci_df, on="LIA_ID")

# Remove rows with zero values and NaNs
filtered_df_flow_velocity = merged_df_flow_velocity[(merged_df_flow_velocity["MEAN_x"] != 0) & (~pd.isna(merged_df_flow_velocity["u_sa"]))]

# Drop rows with NaN values in either column
filtered_df_flow_velocity = filtered_df_flow_velocity.dropna(subset=["MEAN_x", "u_sa"])
filtered_df_flow_velocity["diff_vel"] = filtered_df_flow_velocity["MEAN_x"]-filtered_df_flow_velocity["u_sa"]
filtered_df_flow_velocity["AREA_km2"] = filtered_df_flow_velocity["AREA_x"]/1000000

flow_diff = filtered_df_flow_velocity["diff_vel"]
Area_km2 = filtered_df_flow_velocity["AREA_km2"]
flow_obs = filtered_df_flow_velocity["MEAN_x"]
flow_param = filtered_df_flow_velocity["u_sa"]

slope, intercept, r_value, p_value, std_err = linregress(flow_obs, flow_param)
trendline = slope * flow_obs + intercept

# Scatter plot with trendline
plt.figure()
plt.scatter(flow_obs, flow_param, marker="o", edgecolors="black", facecolors="none", alpha=0.5, s = 10, label="Glaciers")
plt.plot(flow_obs, trendline, color='red', label=f'Linear Trendline: y = {slope:.2f}x + {intercept:.2f}\nR² = {r_value**2:.4f}')
plt.xlabel("Flow velocity 2015-2021 from Rabatel et al. 2023 (m a$^{-1}$)", fontsize = 12)
plt.ylabel("Flow velocity from parameter for 2015 (m a$^{-1}$)", fontsize = 12)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.ylim(0, 100)
plt.xlim(0, 100)
plt.legend()
plt.show()


slope, intercept, r_value, p_value, std_err = linregress(Area_km2,flow_diff)
trendline = slope * flow_obs + intercept
plt.figure()
plt.scatter(Area_km2,flow_diff, marker="o", edgecolors="black", facecolors="none",s = 10, alpha=0.5, label="Glaciers")
#plt.plot(Area_km2, trendline, color='red', label=f'Linear Trendline: y = {slope:.2f}x + {intercept:.2f}\nR² = {r_value**2:.4f}')
plt.xlabel("Area 2015 (km$^{2}$)", fontsize = 12)
plt.ylabel("Difference between methods (m a$^{-1}$)", fontsize = 12)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.xscale("log")  # Set x-axis to logarithmic scale
plt.xlim(0, 100)
plt.ylim(-200, 100)
plt.legend()
plt.show()

#Part IV: export tables
# %% export result dataframes
RGI_df.to_excel("RGI_results_output_new.xlsx", index=False)
cci_df.to_excel("cci_results_output_new.xlsx", index=False)
LIA_df.to_excel("LIA_results_output_new.xlsx", index=False)
SxS_df.to_excel("SxS_results_output_new.xlsx", index=False)
Combined_df.to_excel("Combined_results_output.xlsx", index=False)
combined_df_methods.to_excel("Combined_results_methods_output.xlsx", index=False)

# Export 'region_summary' to an Excel file
region_summary_LIA.to_excel("Region_Summary_LIA.xlsx", index=False)
region_summary_cci = (
    cci_df.groupby("region_id")
    .agg({"V_km3": "sum", "AREA": "sum", "MIN": "mean"})
    .reset_index()
)

# Rename the columns for clarity
region_summary_cci.rename(
    columns={
        "V_km3_cci": "Total_V_km3_cci",
        "AREA_LIA": "Total_AREA_cci",
        "MIN_cci": "Average_MIN_cci",
    },
    inplace=True,
)

# Export 'region_summary' to an Excel file
region_summary_cci.to_excel("Region_Summary_cci.xlsx", index=False)
