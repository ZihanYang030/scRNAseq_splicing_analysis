from pathlib import Path
import os
import re
import time
import pandas as pd
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeoutError

INPUT_FILE = os.environ.get("SRR_LIST", "GSM_OR_RUN_SRR_Acc_List.txt")
OUT_CSV = "original_fastq_s3_names.csv"
# added xlsx file output
OUT_XLSX = "original_fastq_s3_names.xlsx"
FAILED_FILE = "failed_srr.txt"

HEADLESS = True
SLEEP_SECONDS = 1.0 # might need to be adjusted
MAX_SRRS = 5


def read_srr_list(path):
    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip().startswith("SRR")]


def extract_s3_paths_from_page_text(text):
    pattern = r"s3://[^\s\"'<>]+?\.fastq(?:\.gz)?(?:\.\d+)?"
    paths = re.findall(pattern, text)
    cleaned = [p.rstrip(".,;)") for p in paths]
    seen = set()
    unique = []
    for p in cleaned:
        if p not in seen:
            unique.append(p)
            seen.add(p)
    return unique


def infer_read_pair(file_name):
    if re.search(r"(^|[_\-.])R1([_\-.]|$)", file_name):
        return "R1"
    if re.search(r"(^|[_\-.])R2([_\-.]|$)", file_name):
        return "R2"
    return ""


def infer_compression(file_name):
    if ".fastq.gz" in file_name:
        return "gzip"
    if ".fastq" in file_name:
        return "uncompressed_fastq"
    return ""


def scrape_one_srr(page, srr):
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/index.html?run={srr}"
    page.goto(url, wait_until="networkidle", timeout=60000)
    page.wait_for_timeout(3000)

    clicked = False
    for locator_expr in ["text=Data Access", "role=tab[name='Data Access']", "role=button[name='Data Access']"]:
        try:
            if locator_expr.startswith("text="):
                loc = page.locator(locator_expr).first
            elif locator_expr.startswith("role=tab"):
                loc = page.get_by_role("tab", name="Data Access")
            else:
                loc = page.get_by_role("button", name="Data Access")
            if loc.count() > 0:
                loc.click(timeout=5000)
                page.wait_for_timeout(3000)
                clicked = True
                break
        except Exception:
            pass

    page.mouse.wheel(0, 2000)
    page.wait_for_timeout(1000)
    text = page.locator("body").inner_text(timeout=30000)
    s3_paths = extract_s3_paths_from_page_text(text)

    rows = []
    for p in s3_paths:
        file_name = p.split("/")[-1]
        rows.append({
            "SRR": srr, "sra_run_browser_url": url,
            "s3_path": p, "file_name": file_name,
            "read_pair": infer_read_pair(file_name),
            "compression": infer_compression(file_name),
            "data_access_clicked": clicked,
        })
    return rows


def save_outputs(rows, failed):
    columns = ["SRR", "sra_run_browser_url", "s3_path", "file_name", "read_pair", "compression", "data_access_clicked"]
    df = pd.DataFrame(rows, columns=columns).drop_duplicates()
    #if len(df) = 0:
        #break
    if len(df) > 0:
        df = df.sort_values(["SRR", "read_pair", "file_name"])
    df.to_csv(OUT_CSV, index=False)
    df.to_excel(OUT_XLSX, index=False)
    Path(FAILED_FILE).write_text("\n".join(failed) + ("\n" if failed else ""))
    print(f"Rows saved: {len(df)}; failed SRRs: {len(failed)}")


def main():
    srrs = read_srr_list(INPUT_FILE)
    if MAX_SRRS is not None:
        srrs = srrs[:MAX_SRRS]

    all_rows = []
    failed = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=HEADLESS)
        context = browser.new_context(viewport={"width": 1400, "height": 1000})
        page = context.new_page()

        for i, srr in enumerate(srrs, start=1):
            print(f"[{i}/{len(srrs)}] Scraping {srr}")
            try:
                rows = scrape_one_srr(page, srr)
                if rows:
                    all_rows.extend(rows)
                else:
                    failed.append(srr)
            except PlaywrightTimeoutError:
                failed.append(srr)
            except Exception as e:
                print(f"Error for {srr}: {e}")
                failed.append(srr)

            if i % 20 == 0:
                save_outputs(all_rows, failed)
            time.sleep(SLEEP_SECONDS)

        browser.close()

    save_outputs(all_rows, failed)


if __name__ == "__main__":
    main()
