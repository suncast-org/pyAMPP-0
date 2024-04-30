#!/usr/bin/env python3

import os
from sunpy.net import jsoc, fido_factory, Fido, attrs as a
import drms

def download_magnetograms(start_date, end_date, download_folder, jsoc_email, streams=5, print_only=False):
    notifier = a.jsoc.Notify(jsoc_email)
    series = a.jsoc.Series.hmi_b_720s
    segments = a.jsoc.Segment('field') & a.jsoc.Segment('inclination') & a.jsoc.Segment('azimuth')\
        & a.jsoc.Segment('disambig')
    
    series_m = a.jsoc.Series.hmi_m_720s
    segments_m = a.jsoc.Segment("magnetogram")
    
    series_ic = a.jsoc.Series.hmi_ic_nolimbdark_720s
    segments_ic = a.jsoc.Segment("continuum")

    c = drms.Client()
    si = c.info(series.value).segments.index.values
    search_times = a.Time(start_date, end_date)

    search    = Fido.search(search_times, series, notifier, segments)
    search_m  = Fido.search(search_times, series_m, notifier, segments_m)
    search_ic = Fido.search(search_times, series_ic, notifier, segments_ic)

    print(search)
    print(search_m)
    print(search_ic)

    if print_only:
        return

    if not os.path.exists(download_folder):
        os.makedirs(download_folder)

    fetch_params = dict(path=download_folder, max_conn=streams)

    results    = Fido.fetch(search,    **fetch_params)
    results_m  = Fido.fetch(search_m,  **fetch_params)
    results_ic = Fido.fetch(search_ic, **fetch_params)
    results_combined = results + results_m + results_ic

    return results_combined