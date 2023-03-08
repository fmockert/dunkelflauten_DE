"""Combines solar, onshore-wind and offshore-wind capacity factors according to weights."""

__author__      = "Fabian Mockert, Fabian Neumann"
__copyright__   = "Copyright 2019-2020, Fabian Mockert, Fabian Neumann"

import pandas as pd
import numpy as np

import logging
logging.basicConfig(level="INFO")

if __name__ == "__main__":

    config = snakemake.config["technology"]
    country = "DE"#snakemake.wildcards.country

    onwind_w = config['onwind']['weight']
    offwind_w = config['offwind']['weight']
    solar_w = config['solar']['weight']

    onwind_cf = pd.read_csv(snakemake.input.onwind, index_col=0).squeeze()
    offwind_cf = pd.read_csv(snakemake.input.offwind, index_col=0).squeeze()
    solar_cf = pd.read_csv(snakemake.input.solar, index_col=0).squeeze()

    carriers = {"onwind": {"cf": onwind_cf, "w": onwind_w},
                "offwind": {"cf": offwind_cf, "w": offwind_w},
                "solar": {"cf": solar_cf, "w": solar_w}}

    year_start, year_end = snakemake.config["years"]
    idx = pd.date_range(str(year_start), str(year_end+1), freq='H', closed='left')
    weighting_factor = pd.Series(0., index=idx, name=country, dtype="float")
    weighting_sum = 0
    for carrier in carriers.keys():
        if not carriers[carrier]["cf"].empty:
            #weighting_factor += carriers[carrier]["cf"]*carriers[carrier]["w"]
            weighting_factor = pd.Series(weighting_factor.values+(carriers[carrier]["cf"]*carriers[carrier]["w"]).values, index=idx, name=country, dtype="float")
            weighting_sum += carriers[carrier]["w"]
            #print(weighting_factor)
        else:
            logging.info(f"The combined capacity factor for {country} excludes {carrier}.")

    combined = weighting_factor / weighting_sum
    print(combined)
    if combined.isnull().all():
        logging.warning(f"The combined capacity factor is for all values NaN and therefore not generated for {country}. This might raise errors in further processes.")
    else:
        combined.to_csv(snakemake.output[0], header=True)
