import urllib.request

def fetch_pdb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    filename = f"{pdb_id}.pdb"
    print(f"Fetching {pdb_id}...")
    urllib.request.urlretrieve(url, filename)
    print(f"Saved {filename}")

if __name__ == "__main__":
    targets = ["6F86", "4KMU"]
    for target in targets:
        fetch_pdb(target)