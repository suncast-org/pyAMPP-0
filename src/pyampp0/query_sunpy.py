#!/usr/bin/env python3

import os
from sunpy.net import jsoc, fido_factory, Fido, attrs as a
from astropy.time import Time
import numpy as np
import datetime
import drms

def download_closest_magnetograms(target_datetime, download_folder, jsoc_email, streams=5):
    """Downloads the closest magnetogram to a given datetime.

    Args:
        target_datetime (datetime): The target date and time for which to find the closest 
                                     magnetogram.
        download_folder (str): Path to the folder where downloaded files will be saved.
        jsoc_email (str): Your JSOC email address for authentication.
        streams (int, optional): Number of concurrent downloads (defaults to 5).

    Returns:
        sunpy.net.Fido.FileCollection or None: A collection of the downloaded files
                                             if successful; None otherwise.
    """

    notifier = a.jsoc.Notify(jsoc_email)
    series = a.jsoc.Series.hmi_b_720s
    segments = a.jsoc.Segment('field') & a.jsoc.Segment('inclination') & a.jsoc.Segment('azimuth') \
        & a.jsoc.Segment('disambig')

    series_m = a.jsoc.Series.hmi_m_720s
    segments_m = a.jsoc.Segment("magnetogram")
    
    series_ic = a.jsoc.Series.hmi_ic_nolimbdark_720s
    segments_ic = a.jsoc.Segment("continuum")

    c = drms.Client()
    si = c.info(series.value).segments.index.values

    # choosing the closest time
    start_date = target_datetime - datetime.timedelta(minutes=6)
    end_date   = target_datetime + datetime.timedelta(minutes=6)
    search_times = a.Time(start_date, end_date)
    search    = Fido.search(search_times, series, notifier, segments)

    target_atime = Time(target_datetime)
    times_list = [Time(datetime.datetime.strptime(x, "%Y.%m.%d_%H:%M:%S_TAI")) for x in search["jsoc"]["T_REC"]]
    times_deltas = np.array([np.abs((x-target_atime).sec) for x in times_list])
    closest_pos = np.argmin(times_deltas)
    closest_time = times_list[closest_pos]
    closest_timestr = search["jsoc"]["T_REC"][closest_pos]
    
    # querying the right magnetograms

    search_times = a.Time(closest_time, closest_time) # Set search to single closest time

    search = Fido.search(search_times, series, notifier, segments)
    search_m = Fido.search(search_times, series_m, notifier, segments_m)
    search_ic = Fido.search(search_times, series_ic, notifier, segments_ic)

    #print(search)
    #print(search_m)
    #print(search_ic)

    #if print_only:
    #    return
    print(f"downloading for time {closest_time.iso}")

    if not os.path.exists(download_folder):
        os.makedirs(download_folder)

    fetch_params = dict(path=download_folder, max_conn=streams)

    results = Fido.fetch(search, **fetch_params)
    results_m = Fido.fetch(search_m, **fetch_params)
    results_ic = Fido.fetch(search_ic, **fetch_params)
    return results + results_m + results_ic 
