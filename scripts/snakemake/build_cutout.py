"""Downloads reanalysis weather data e.g. from ERA-5."""

__author__      = "Fabian Mockert, Fabian Neumann"
__copyright__   = "Copyright 2023, Fabian Mockert, Fabian Neumann"

import os
import atlite
import sys

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level="INFO")

if __name__ == "__main__":

    cutout_params = snakemake.config["cutout"]

    year_start, year_end = cutout_params.pop("years")
    logging.info(f"Download weather data from {year_start} to {year_end}.")

    x0, x1 = cutout_params["xs"]
    y0, y1 = cutout_params["ys"]
    logging.info(f"Rectangle is x0={x0} x1={x1} y0={y0} y1={y1}")

    for i in range(year_start, year_end+1):

        logging.info(f"Start downloading year {i}.")



        cutout = atlite.Cutout(
                    module='era5', path=os.path.dirname(snakemake.output[0]+'cutouts_germany_%s.nc'%(i)),
                    x=slice(x0,x1), y=slice(y0,y1),
                    years=slice(i,i))

        cutout.prepare()
