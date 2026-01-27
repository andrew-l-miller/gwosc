import requests
#from gwpy.timeseries import TimeSeries


def fetch_run_gps_times(run):
    "Fetch gwosc archive and return the (start, end) GPS time tuple of the run."
    response = requests.get("https://gwosc.org/archive/all/json/").json()
    runs = response["runs"]
    run_info = runs.get(run)
    if run_info is None:
        raise ValueError(f"Could not find run {run}. Available runs: {runs.keys()}")
    return run_info["GPSstart"], run_info["GPSend"]


def fetch_strain_list(run, detector, gps_start=None, gps_end=None):
    "Return the list of strain file info for `run` and `detector`."
    if gps_start is None or gps_end is None:
        start, end = fetch_run_gps_times(run)
        gps_start = gps_start or start
        gps_end = gps_end or end

    # Get the strain list
    fetch_url = (
        f"https://gwosc.org/archive/links/"
        f"{run}/{detector}/{gps_start}/{gps_end}/json/"
    )
    response = requests.get(fetch_url)
    response.raise_for_status()
    return response.json()["strain"]


def download_strain_file(download_url,outdir='./'):
    "Download the strain file on the given url and save to disk."
    # In the next line I parse the file name from the download url.
    # Ideally, the file name should be grabbed from the
    # Content-Disposition response header.
    filename = download_url.split("/")[-1]
    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with open(outdir+filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return filename




import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Download GWOSC strain data")
    parser.add_argument("channel", nargs="?", default="O3b_4KHZ_R1", help="Data channel (default: O3b_4KHZ_R1)")
    parser.add_argument("detector", nargs="?", default="H1", help="Detector (H1 or L1, default: H1)")
    parser.add_argument("save_dir", nargs="?", default="./data/", help="Directory to save files (default: ./data)")
    parser.add_argument("--gps-start", type=int, default=None,
                        help="GPS start time (default: run start)")
    parser.add_argument("--gps-end", type=int, default=None,
                        help="GPS end time (default: run end)")

    args = parser.parse_args()

    # make sure directory exists
    os.makedirs(args.save_dir, exist_ok=True)

    # fetch strain files
    strain_files = fetch_strain_list(args.channel, args.detector,gps_start=args.gps_start,gps_end=args.gps_end)
    try:
        with open("filesdone.txt", "r") as fp:
            donelist = [f.strip() for f in fp.readlines()]
    except FileNotFoundError:
        donelist = []

    for afile in strain_files:
        if afile["url"] in donelist:
            continue
        if afile["format"] == "gwf":
            print(f"Downloading {afile['url']}")
            fname = download_strain_file(afile["url"], outdir=args.save_dir)
            # tseries = TimeSeries.read(fname, format="hdf5.gwosc")
            with open("filesdone.txt", "a") as fp:
                fp.write(f"{afile['url']}\n")
            # process tseries here


if __name__ == "__main__":
    main()  
