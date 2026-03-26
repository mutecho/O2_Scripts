#!/usr/bin/env python3
import requests
import json
import sys

def check_run(run_number):
    url = f"http://alice-ccdb.cern.ch/browse/RCT/Info/RunInformation/{run_number}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return run_number, "HTTP Error", r.status_code
        data = r.json()
        if "objects" not in data or not data["objects"]:
            return run_number, "No data", None
        obj = data["objects"][0]
        size = obj.get("size", 0)
        ctype = obj.get("contentType", "unknown")
        if size <= 1:
            return run_number, "Invalid (size<=1)", f"{size}, {ctype}"
        return run_number, "OK", f"size={size}, type={ctype}"
    except Exception as e:
        return run_number, "Error", str(e)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python check_ccdb_runs.py <start_run> <end_run>")
        sys.exit(1)

    start_run = int(sys.argv[1])
    end_run = int(sys.argv[2])

    print(f"检查 CCDB Run 信息区间: {start_run}–{end_run}")
    print("RunNumber\tStatus\t\tDetails")

    for run in range(start_run, end_run + 1):
        run_no, status, detail = check_run(run)
        print(f"{run_no}\t{status:15}\t{detail}")