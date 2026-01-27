import sys
from functions_to_download_data import *

if len(sys.argv) < 1:
    print("Usage: python3 function.py <save_dir>")
    sys.exit(1)

save_dir = sys.argv[1]


strain_files = fetch_strain_list("O3b_4KHZ_R1", "H1", 1264525312, 1264550417)
print(f"Found {len(strain_files)} files")
#print(strain_files[0:5])

for a_file in strain_files:
    if a_file["GPSstart"] == 1264525312 and a_file["format"] == "gwf":
        print(f"Downloading {a_file['url']}")
        fname = download_strain_file(a_file["url"],save_dir)
