import requests


def fetch_strain_list(run, detector, gps_start, gps_end):
    "Return the list of strain file info for `run` and `detector`."

    # Get the strain list
    fetch_url = (
        f"https://gwosc.org/archive/links/"
        f"{run}/{detector}/{gps_start}/{gps_end}/json/"
    )
    response = requests.get(fetch_url)
    response.raise_for_status()
    return response.json()["strain"]


def download_strain_file(download_url,savedir='./'):
    "Download the strain file on the given url and save to disk."
    # In the next line I parse the file name from the download url.
    # Ideally, the file name should be grabbed from the
    # Content-Disposition response header.
    filename = download_url.split("/")[-1]
    with requests.get(download_url, stream=True) as r:
        r.raise_for_status()
        with open(savedir+filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return filename
