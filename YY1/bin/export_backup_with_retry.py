#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import subprocess

TOKENS = {
    "YY1-CNVseq汇总.xlsx": "NWgXs0zzShsXYLtPuN5c1M3fn52",
    "YY1-PGTA汇总.xlsx": "UfFEsGG0Zh97KotvvWzcLRsNnSf",
    "YY1-49SD汇总.xlsx": "FbLRsZf6OhGceSthhflcKSVrnIe",
    "YY1-tNGS汇总.xlsx": "YLmusKmiphU8iPtTz6zc5bo1ni9",
    "YY1-mNGS汇总.xlsx": "MjihsiLrnh6UU5tc4cSc4PRSntc",
    "YY1-TB2汇总.xlsx": "N6FRsWZGRhaLskt57hTckWFlnSc"
}

BACKUP_DIR = "/mnt/gpfs1/Dataset/04.project/Bioinfo/99.Sequencer_assessment/01.BK/02.Analysis/04.new_qican/S100/SE75/summary/YY1"
LARK_CLI = "/mnt/gpfs1/Users/yangjinxurong/software/lark-cli/bin/lark-cli"

def export_with_retry(name, token):
    output_path = f"./{name}"
    cmd = f"{LARK_CLI} sheets +workbook-export --spreadsheet-token '{token}' --output-path '{output_path}'"
    
    max_retries = 8
    for attempt in range(1, max_retries + 1):
        print(f"[INFO] Exporting {name} (Attempt {attempt}/{max_retries})...")
        res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0:
            print(f"[SUCCESS] Exported {name} successfully.")
            return True
        else:
            print(f"[WARN] Attempt {attempt} failed:\nStdout: {res.stdout}\nStderr: {res.stderr}")
            if attempt < max_retries:
                print("[INFO] Sleeping 10 seconds before next retry...")
                time.sleep(10)
    print(f"[ERROR] Failed to export {name} after {max_retries} attempts.")
    return False

def main():
    print("=== Starting Workbook Backup with Retry Logic ===")
    os.chdir(BACKUP_DIR)
    
    print("[INFO] Waiting 15 seconds for Feishu backend to stabilize...")
    time.sleep(15)
    
    success_count = 0
    for name, token in TOKENS.items():
        if export_with_retry(name, token):
            success_count += 1
            
    print(f"=== Backup completed: {success_count}/{len(TOKENS)} succeeded ===")
    if success_count < len(TOKENS):
        sys.exit(1)

if __name__ == "__main__":
    main()
