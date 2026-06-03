#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
scrape_kahaku_azami_all_selenium.py

科博「日本のアザミ」DBをブラウザ（Selenium）で開いて、レンダリング後の全記録を保存する。
requests で detail.html がテンプレだけ返る場合の回避用。

出力:
  kahaku_azami_all_records_selenium.csv
  kahaku_azami_all_records_selenium_rawtext.csv
  kahaku_html_cache/no_XXX.html
  kahaku_screenshots/no_XXX.png   # --screenshots を付けた場合

使い方:
  python -m pip install selenium pandas beautifulsoup4
  python scrape_kahaku_azami_all_selenium.py --start 1 --end 220

ブラウザを見ながら確認したい場合:
  python scrape_kahaku_azami_all_selenium.py --start 1 --end 220 --no-headless

スクショも保存:
  python scrape_kahaku_azami_all_selenium.py --start 1 --end 220 --screenshots
"""

import argparse
import re
import time
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup

from selenium import webdriver
from selenium.webdriver.chrome.options import Options as ChromeOptions
from selenium.webdriver.common.by import By


BASE = "https://www.kahaku.go.jp/research/db/botany/azami/"
DETAIL = BASE + "detail.html?no={no}"

FIELDS = [
    "分類", "和名", "別名", "種名", "変種名", "キャッチフレーズ",
    "基準産地", "記載", "分布", "ノート", "撮影データ", "分布図"
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--start", type=int, default=1)
    p.add_argument("--end", type=int, default=250)
    p.add_argument("--out", default="kahaku_azami_all_records_selenium.csv")
    p.add_argument("--sleep-sec", type=float, default=1.0)
    p.add_argument("--wait-sec", type=float, default=2.0)
    p.add_argument("--no-headless", action="store_true")
    p.add_argument("--screenshots", action="store_true")
    p.add_argument("--cache-html", action="store_true", default=True)
    return p.parse_args()


def make_driver(headless=True):
    opts = ChromeOptions()
    if headless:
        opts.add_argument("--headless=new")
    opts.add_argument("--window-size=1400,1800")
    opts.add_argument("--disable-gpu")
    opts.add_argument("--lang=ja-JP")
    opts.add_argument("--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/124 Safari/537.36")

    # Selenium Manager should automatically find/download chromedriver in recent Selenium.
    driver = webdriver.Chrome(options=opts)
    return driver


def clean_text(x):
    if x is None:
        return ""
    return re.sub(r"\s+", " ", str(x)).strip()


def lines_from_text(text):
    raw = [x.strip() for x in text.splitlines()]
    return [x for x in raw if x]


def extract_field_from_lines(lines, field):
    """
    visible text の行列から field の次の値を取る。
    値が複数行の場合は次の field まで結合。
    """
    out = ""
    for i, line in enumerate(lines):
        if line == field:
            vals = []
            for j in range(i + 1, len(lines)):
                if lines[j] in FIELDS:
                    break
                # navigation/footer noise
                if lines[j].startswith("←") or "前ページ" in lines[j] or "トップへ" in lines[j]:
                    break
                vals.append(lines[j])
            out = clean_text(" ".join(vals))
            break
    return out


def extract_table_like_fields(soup):
    """
    th/td や dt/dd のような構造があれば拾う。
    """
    d = {}

    # th/td rows
    for tr in soup.find_all("tr"):
        cells = [clean_text(c.get_text(" ")) for c in tr.find_all(["th", "td"])]
        if len(cells) >= 2:
            key = cells[0]
            val = clean_text(" ".join(cells[1:]))
            if key in FIELDS and val:
                d[key] = val

    # dt/dd
    for dt in soup.find_all("dt"):
        key = clean_text(dt.get_text(" "))
        dd = dt.find_next_sibling("dd")
        if key in FIELDS and dd:
            val = clean_text(dd.get_text(" "))
            if val:
                d[key] = val

    return d


def extract_species(text):
    m = re.search(r"\b(Cirsium\s+[a-z][a-z\-]+)\b", text)
    return m.group(1) if m else ""


def extract_orientation(text):
    # text should include catchphrase/notes
    chunks = re.split(r"[。．\.]", text)
    phrases = []
    for ch in chunks:
        if any(k in ch for k in ["頭花", "花を"]):
            if any(k in ch for k in ["上向", "斜上", "直立", "下向", "点頭", "垂下", "垂れ"]):
                phrases.append(clean_text(ch))

    phrase = " / ".join(phrases[:4])

    target = phrase if phrase else text

    up = bool(re.search(r"(頭花|花)[^。．\.]{0,120}(上向き|上向|斜上|直立)", target))
    down = bool(re.search(r"(頭花|花)[^。．\.]{0,120}(下向き|下向|点頭|垂下|垂れ)", target))

    # Catch specific wording like "上向きないし斜め"
    if "上向きないし斜め" in target or "斜め上向き" in target:
        up = True

    if up and down:
        return "mixed_or_review", phrase
    if up:
        return "upward", phrase
    if down:
        return "nodding", phrase
    return "unknown", phrase


def scrape_one(driver, no, wait_sec=2.0, save_html_dir=None, screenshot_dir=None):
    url = DETAIL.format(no=no)
    driver.get(url)
    time.sleep(wait_sec)

    html = driver.page_source
    body_text = driver.find_element(By.TAG_NAME, "body").text if driver.find_elements(By.TAG_NAME, "body") else ""
    soup = BeautifulSoup(html, "html.parser")

    if save_html_dir is not None:
        save_html_dir.mkdir(exist_ok=True)
        (save_html_dir / f"no_{no:03d}.html").write_text(html, encoding="utf-8", errors="ignore")

    if screenshot_dir is not None:
        screenshot_dir.mkdir(exist_ok=True)
        try:
            driver.save_screenshot(str(screenshot_dir / f"no_{no:03d}.png"))
        except Exception:
            pass

    # collect text from visible body + meta + alt
    extra = []
    for tag in soup.find_all("meta"):
        val = tag.get("content")
        if val:
            extra.append(val)
    for tag in soup.find_all(["img", "a"]):
        for attr in ["alt", "title"]:
            val = tag.get(attr)
            if val:
                extra.append(val)

    all_text = clean_text(body_text + "\n" + "\n".join(extra))
    lines = lines_from_text(body_text)

    field_from_table = extract_table_like_fields(soup)

    rec = {"no": no, "url": url}
    for f in FIELDS:
        rec[f] = field_from_table.get(f, "") or extract_field_from_lines(lines, f)

    # fallback species/japanese name
    species = rec.get("種名", "") or extract_species(all_text)
    rec["species_extracted"] = species

    japanese = rec.get("和名", "")
    if not japanese:
        # try before species
        if species and species in all_text:
            before = all_text[max(0, all_text.find(species) - 80):all_text.find(species)]
            ms = list(re.finditer(r"([ァ-ヶー一-龥・（）\(\)]{2,35}アザミ(?:（新称）)?)", before))
            if ms:
                japanese = ms[-1].group(1)
        if not japanese:
            m = re.search(r"([ァ-ヶー一-龥・（）\(\)]{2,35}アザミ(?:（新称）)?)", all_text)
            if m and "図鑑" not in m.group(1):
                japanese = m.group(1)
    rec["japanese_name_extracted"] = japanese

    orientation_text_source = " ".join([
        rec.get("キャッチフレーズ", ""),
        rec.get("ノート", ""),
        rec.get("記載", ""),
        all_text
    ])
    orientation, phrase = extract_orientation(orientation_text_source)
    rec["head_orientation_binary"] = orientation
    rec["orientation_phrase"] = phrase

    # status
    if not species and not japanese and not phrase:
        rec["record_status"] = "no_data_visible"
    elif not species and not phrase:
        rec["record_status"] = "weak_data"
    else:
        rec["record_status"] = "ok"

    rec["visible_text_preview"] = clean_text(body_text)[:1200]
    rec["all_text_preview"] = all_text[:1600]
    return rec


def main():
    args = parse_args()

    out_path = Path(args.out)
    raw_out = out_path.with_name(out_path.stem + "_rawtext.csv")
    html_dir = Path("kahaku_html_cache") if args.cache_html else None
    screenshot_dir = Path("kahaku_screenshots") if args.screenshots else None

    driver = make_driver(headless=not args.no_headless)

    rows = []
    try:
        for no in range(args.start, args.end + 1):
            try:
                rec = scrape_one(
                    driver,
                    no,
                    wait_sec=args.wait_sec,
                    save_html_dir=html_dir,
                    screenshot_dir=screenshot_dir
                )
            except Exception as e:
                rec = {
                    "no": no,
                    "url": DETAIL.format(no=no),
                    "record_status": "error",
                    "error": repr(e)
                }

            rows.append(rec)
            print(
                f"no={no:03d} status={rec.get('record_status')} "
                f"species={rec.get('species_extracted','')} "
                f"name={rec.get('japanese_name_extracted','')} "
                f"orientation={rec.get('head_orientation_binary','')}"
            )

            # incremental save
            df = pd.DataFrame(rows)
            df.to_csv(out_path, index=False, encoding="utf-8-sig")
            time.sleep(args.sleep_sec)

    except KeyboardInterrupt:
        print("\nInterrupted. Saved partial output.")
    finally:
        try:
            driver.quit()
        except Exception:
            pass

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False, encoding="utf-8-sig")

    raw_cols = ["no", "url", "record_status", "visible_text_preview", "all_text_preview"]
    raw_cols = [c for c in raw_cols if c in df.columns]
    df[raw_cols].to_csv(raw_out, index=False, encoding="utf-8-sig")

    summary = (
        df.groupby(["record_status", "head_orientation_binary"], dropna=False)
          .size()
          .reset_index(name="n")
          .sort_values(["record_status", "head_orientation_binary"])
    )
    summary_out = out_path.with_name(out_path.stem + "_summary.csv")
    summary.to_csv(summary_out, index=False, encoding="utf-8-sig")

    print("\nSaved:")
    print(" ", out_path)
    print(" ", raw_out)
    print(" ", summary_out)
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
