from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "_site"
EDITION = (ROOT / "COURSE_EDITION").read_text(encoding="utf-8").strip()

if not SITE.exists():
    raise SystemExit("_site does not exist; run quarto render first")


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self.hrefs.append(href)


broken: list[tuple[Path, str]] = []
html_files = list(SITE.rglob("*.html"))

for page in html_files:
    text = page.read_text(encoding="utf-8")
    parser = LinkParser()
    parser.feed(text)

    if "<p>[" in text or "<p> [" in text:
        raise SystemExit(f"Possible malformed display math in {page.relative_to(SITE)}")

    for href in parser.hrefs:
        if href.startswith(("#", "http://", "https://", "mailto:", "javascript:")):
            continue
        target_text = unquote(href.split("#", 1)[0].split("?", 1)[0])
        if not target_text:
            continue
        target = (page.parent / target_text).resolve()
        if target.is_dir():
            target = target / "index.html"
        if not target.exists():
            broken.append((page.relative_to(SITE), href))

if broken:
    details = "\n".join(f"{page}: {href}" for page, href in broken)
    raise SystemExit(f"Broken rendered links:\n{details}")

required_pages = [
    SITE / "index.html",
    SITE / "modules/01-count-matrix/index.html",
    SITE / "modules/02-quality-control/index.html",
    SITE / "modules/03-normalization/index.html",
    SITE / "labs/index.html",
]
missing_pages = [str(path.relative_to(SITE)) for path in required_pages if not path.exists()]
if missing_pages:
    raise SystemExit(f"Missing rendered pages: {missing_pages}")

if EDITION == "student":
    public_forbidden = [
        SITE / "instructor",
        SITE / "solutions",
        SITE / "rubrics",
        SITE / "teaching-notes",
        SITE / "expected-results",
    ]
    present_forbidden = [
        str(path.relative_to(SITE)) for path in public_forbidden if path.exists()
    ]
    if present_forbidden:
        raise SystemExit(
            f"Instructor-only content rendered publicly: {present_forbidden}"
        )

print(
    f"Rendered-site validation passed: {len(html_files)} HTML pages and all "
    f"local links are valid for the {EDITION} edition."
)
