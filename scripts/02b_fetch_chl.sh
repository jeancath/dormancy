#!/bin/sh
# ESA OC-CCI v6.0 monthly chlorophyll-a, 1998-2022, subsampled every 6th 1/24-deg pixel (~0.25 deg)
# from NOAA PIFSC ERDDAP. One NetCDF per year.
for y in $(seq 1998 2022); do
  f=chl_$y.nc
  [ -s $f ] && continue
  curl -s -m 1800 --retry 3 -o $f "https://oceanwatch.pifsc.noaa.gov/erddap/griddap/esa-cci-chla-monthly-v6-0.nc?chlor_a%5B(${y}-01-01T00:00:00Z):1:(${y}-12-01T00:00:00Z)%5D%5B(89.97916666666667):6:(-89.97916666666666)%5D%5B(0.020833333333314386):6:(359.97916666666663)%5D"
  echo "$y $(stat -c %s $f)"
done
echo CHLDONE
