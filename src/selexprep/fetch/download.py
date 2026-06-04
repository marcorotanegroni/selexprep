"""Download FASTQ data for SELEX BioProjects.

**v0.1 backend order (MIT-compatible by default):**

1. ``ENA-direct`` via ``requests`` — default; pulls FASTQ + published MD5
   from ENA's ``filereport`` API with Range-resume and bit-level integrity
   verification. No external tools required.
2. ``kingfisher`` (optional subprocess) — used only if installed
   separately and runtime-detected. **GPL-3.0**, NOT bundled with the
   MIT-licensed selexprep package; users install it themselves if they
   want SRA Toolkit / AWS / GCP / Aspera download paths.
3. ``sra-toolkit`` (optional, via ``prefetch`` + ``fasterq-dump``) — used
   only if installed separately.

This reordering (vs the original selex_corpus default of kingfisher-first)
is a licensing requirement: with kingfisher absent, the package must still
work on a clean pip install. Documented behaviour: a default ``pip install
selexprep`` produces downloads via ENA-direct only.

External processed-data sources (Zenodo / Figshare) bypass this backend
chain and use ``requests`` directly with optional MD5 verification.
"""

from __future__ import annotations

import hashlib
import json
import logging
import shutil
import subprocess
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

import requests

from selexprep._common import iter_srr_files

logger = logging.getLogger(__name__)

# Public type alias for the dispatcher's `backend` parameter.
DownloadBackend = Literal["auto", "ena", "kingfisher", "sra"]


# ---------------------------------------------------------------------------
# Filesystem layout helpers
# ---------------------------------------------------------------------------


def safe_dir_name(bioproject_id: str) -> str:
    """Sanitise bioproject_id for use as directory name (colons → underscores)."""
    return bioproject_id.replace(":", "_")


def round_dir(
    raw_root: Path,
    bioproject_id: str,
    round_number: int | None,
    target_hint: str | None = None,
) -> Path:
    """Construct the canonical per-round output directory.

    Layout:
        raw_root/<safe_bp_id>/round_NN/                (single-target)
        raw_root/<safe_bp_id>/<target>/round_NN/       (multi-target)
        raw_root/<safe_bp_id>/round_unknown/           (round not yet assigned)
    """
    base = raw_root / safe_dir_name(bioproject_id)
    if target_hint:
        base = base / target_hint
    if round_number is None:
        return base / "round_unknown"
    return base / f"round_{round_number:02d}"


def srr_present(raw_root: Path, bp_id: str, srr: str) -> bool:
    """True iff at least one canonical FASTQ file for `srr` exists under `raw_root`."""
    bp_raw_dir = raw_root / safe_dir_name(bp_id)
    return bool(iter_srr_files(bp_raw_dir, srr))


def needs_manual_review(round_rec: dict) -> bool:
    """Parse a `needs_manual_review` CSV cell (string truthy → bool)."""
    raw = (round_rec.get("needs_manual_review") or "").strip().lower()
    return raw in ("true", "1", "yes", "y")


# ---------------------------------------------------------------------------
# Download integrity checks
# ---------------------------------------------------------------------------

_MIN_FASTQ_GZ_BYTES = 1024


def validate_fastq_gz(path: Path, min_size: int = _MIN_FASTQ_GZ_BYTES) -> bool:
    """Return True iff `path` exists, is at least `min_size` bytes, and passes ``gzip -t``.

    Generous 600s timeout: a 7-8 GB gzipped FASTQ takes 2-5 min to decompress
    single-threaded; a tighter timeout would false-reject canonical files.
    """
    if not path.exists():
        return False
    try:
        if path.stat().st_size < min_size:
            return False
    except OSError:
        return False
    try:
        result = subprocess.run(
            ["gzip", "-t", str(path)],
            capture_output=True,
            timeout=600,
        )
        if result.returncode != 0:
            stderr_tail = (result.stderr or b"").decode("utf-8", "replace").strip()[-200:]
            logger.warning(
                "    gzip -t failed for %s (rc=%d): %s",
                path.name,
                result.returncode,
                stderr_tail,
            )
            return False
        return True
    except subprocess.TimeoutExpired:
        logger.warning("    gzip -t timed out (>600s) on %s — treating as invalid", path.name)
        return False
    except (subprocess.SubprocessError, OSError) as e:
        logger.warning("    gzip -t error on %s: %s", path.name, e)
        return False


def md5_file(path: Path, chunk: int = 1 << 20) -> str:
    """Stream-MD5 a file in 1 MiB chunks (memory-bounded for multi-GB FASTQs)."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Range-resumable streaming downloader
# ---------------------------------------------------------------------------


def stream_download(
    url: str,
    dest: Path,
    expected_md5: str | None = None,
    timeout: int = 300,
    max_attempts: int = 3,
) -> bool:
    """Stream `url` to `dest` with Range-based resume and optional MD5 verify.

    Idempotent: re-running is cheap when the file is already complete.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    total_size = 0
    try:
        head = requests.head(url, timeout=30, allow_redirects=True)
        if head.ok:
            total_size = int(head.headers.get("Content-Length", 0))
    except requests.RequestException:
        pass

    for attempt in range(1, max_attempts + 1):
        existing = dest.stat().st_size if dest.exists() else 0

        if total_size and existing == total_size:
            logger.info("    [complete] %s (%d bytes)", dest.name, total_size)
            break

        headers: dict[str, str] = {}
        mode = "wb"
        if existing and total_size and existing < total_size:
            headers["Range"] = f"bytes={existing}-"
            mode = "ab"
            logger.info("    resuming %s at byte %d/%d", dest.name, existing, total_size)

        try:
            resp = requests.get(url, headers=headers, stream=True, timeout=timeout)
            resp.raise_for_status()
            with open(dest, mode) as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    if chunk:
                        fh.write(chunk)
        except requests.RequestException as e:
            logger.warning(
                "    [attempt %d/%d] download error for %s: %s",
                attempt,
                max_attempts,
                dest.name,
                e,
            )
            time.sleep(2**attempt)
            continue

        got = dest.stat().st_size if dest.exists() else 0
        if total_size and got != total_size:
            logger.warning(
                "    [attempt %d/%d] size mismatch for %s: got %d, expected %d",
                attempt,
                max_attempts,
                dest.name,
                got,
                total_size,
            )
            time.sleep(2**attempt)
            continue
        break
    else:
        logger.error("    exhausted %d attempts downloading %s", max_attempts, dest.name)
        return False

    if expected_md5:
        got = md5_file(dest)
        if got.lower() != expected_md5.lower():
            logger.error(
                "    md5 mismatch for %s (expected %s, got %s) — deleting",
                dest.name,
                expected_md5,
                got,
            )
            dest.unlink(missing_ok=True)
            return False
        logger.info("    [md5 ok] %s", dest.name)

    return True


# ---------------------------------------------------------------------------
# Backend probes
# ---------------------------------------------------------------------------


def kingfisher_available() -> bool:
    """Probe for the optional ``kingfisher`` subprocess backend (GPL-3.0)."""
    return shutil.which("kingfisher") is not None


def sratoolkit_available() -> bool:
    """Probe for the optional ``prefetch`` + ``fasterq-dump`` subprocess backend."""
    return shutil.which("prefetch") is not None and shutil.which("fasterq-dump") is not None


def _kingfisher_download_methods() -> list[str]:
    methods: list[str] = []
    if shutil.which("ascp"):
        methods.append("ena-ascp")
    if shutil.which("curl"):
        methods.append("ena-ftp")
    if sratoolkit_available():
        methods.append("prefetch")
    if shutil.which("aria2c"):
        methods.append("aws-http")
    if shutil.which("aws"):
        methods.append("aws-cp")
    if shutil.which("gsutil"):
        methods.append("gcp-cp")
    return methods


def _find_sra_cache(output_dir: Path, srr: str) -> Path | None:
    candidates = [
        output_dir / f"{srr}.sra",
        output_dir / srr / f"{srr}.sra",
        *sorted(output_dir.rglob(f"{srr}.sra")),
    ]
    return next((p for p in candidates if p.exists()), None)


def _skip_if_complete_else_clean(srr: str, output_dir: Path) -> bool:
    """Shared backend preamble: make `output_dir`, then decide skip-vs-proceed.

    Returns True iff valid FASTQs for `srr` are already present (caller should
    return True early). Otherwise deletes any partial/corrupt leftovers and
    returns False so the backend proceeds with a clean download.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    fastq = iter_srr_files(output_dir, srr)
    if fastq and all(validate_fastq_gz(p) for p in fastq):
        logger.info("  [skip] %s already present in %s", srr, output_dir)
        return True
    if fastq:
        logger.warning("  %s has partial/corrupt FASTQs in %s — re-downloading", srr, output_dir)
        for p in fastq:
            p.unlink(missing_ok=True)
    return False


# ---------------------------------------------------------------------------
# Backend: ENA-direct (DEFAULT, no external deps)
# ---------------------------------------------------------------------------


def download_srr_ena_direct(srr: str, output_dir: Path, dry_run: bool = False) -> bool:
    """Download an SRR via ENA's filereport API + Range-resumable streaming.

    Default v0.1 backend: no external tools, MIT-compatible. Verifies bytes
    against ENA's published MD5; skips redundant ``gzip -t`` when MD5
    matches.
    """
    if _skip_if_complete_else_clean(srr, output_dir):
        return True

    if dry_run:
        return True

    api_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={srr}&result=read_run"
        "&fields=fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
    )
    try:
        resp = requests.get(api_url, timeout=60)
        resp.raise_for_status()
    except requests.RequestException as e:
        logger.error("  ena-direct: filereport API failed for %s: %s", srr, e)
        return False

    lines = [ln for ln in resp.text.strip().splitlines() if ln]
    if len(lines) < 2:
        logger.error("  ena-direct: filereport returned no rows for %s", srr)
        return False

    header = lines[0].split("\t")
    try:
        ftp_idx = header.index("fastq_ftp")
        md5_idx = header.index("fastq_md5")
    except ValueError:
        logger.error("  ena-direct: unexpected filereport header for %s: %s", srr, header)
        return False

    row = lines[1].split("\t")
    ftp_field = row[ftp_idx] if ftp_idx < len(row) else ""
    md5_field = row[md5_idx] if md5_idx < len(row) else ""
    if not ftp_field:
        logger.error("  ena-direct: filereport returned no fastq_ftp for %s", srr)
        return False

    urls = [u for u in ftp_field.split(";") if u]
    md5s_raw: list[str | None] = list(md5_field.split(";")) if md5_field else []
    if len(md5s_raw) != len(urls):
        md5s_raw = [None] * len(urls)
    md5s: list[str | None] = [m if m else None for m in md5s_raw]

    logger.info("  ena-direct %s → %d file(s) from ENA filereport", srr, len(urls))

    produced: list[Path] = []
    ok_all = True
    for url_path, expected_md5 in zip(urls, md5s, strict=True):
        url = (
            url_path
            if url_path.startswith(("http://", "https://", "ftp://"))
            else f"https://{url_path}"
        )
        dest = output_dir / Path(url_path).name
        if not stream_download(url, dest, expected_md5=expected_md5, timeout=1800, max_attempts=5):
            ok_all = False
            break
        produced.append(dest)

    if not ok_all:
        logger.warning("  ena-direct: cleaning up partial files for %s", srr)
        for p in produced:
            p.unlink(missing_ok=True)
        return False

    md5_verified_all = bool(md5s) and all(m for m in md5s)
    if md5_verified_all:
        logger.info("  ena-direct: all %d file(s) MD5-verified, skipping gzip -t", len(produced))
        return True

    bad = [fq for fq in produced if not validate_fastq_gz(fq)]
    if bad:
        logger.warning(
            "  ena-direct: %d invalid FASTQ(s) for %s after download (%s); deleting all",
            len(bad),
            srr,
            ", ".join(p.name for p in bad),
        )
        for fq in produced:
            fq.unlink(missing_ok=True)
        return False
    return True


# ---------------------------------------------------------------------------
# Backend: kingfisher (optional, GPL-3.0)
# ---------------------------------------------------------------------------


def download_srr_kingfisher(srr: str, output_dir: Path, dry_run: bool = False) -> bool:
    """Download an SRR via the ``kingfisher`` subprocess backend (GPL-3.0).

    Only used when the binary is installed separately. Tries multiple
    download methods (ENA/AWS/GCP/Aspera) depending on which auxiliary
    tools are present.
    """
    if _skip_if_complete_else_clean(srr, output_dir):
        return True

    methods = _kingfisher_download_methods()
    if not methods:
        logger.warning("  kingfisher available but no supported download methods detected")
        return False

    cmd = [
        "kingfisher",
        "get",
        "-r",
        srr,
        "--download-methods",
        *methods,
        "--output-directory",
        str(output_dir),
        "--output-format-possibilities",
        "fastq.gz",
        "--download-threads",
        "4",
        "--check-md5sums",
    ]
    logger.info("  kingfisher %s → %s (methods: %s)", srr, output_dir, ",".join(methods))
    if dry_run:
        return True

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            logger.warning("  kingfisher failed for %s: %s", srr, result.stderr[-500:])
            return False
    except subprocess.TimeoutExpired:
        logger.error("  kingfisher timeout for %s", srr)
        return False

    produced = iter_srr_files(output_dir, srr)
    if not produced:
        logger.warning("  kingfisher claimed success for %s but no FASTQ found", srr)
        return False
    bad = [fq for fq in produced if not validate_fastq_gz(fq)]
    if bad:
        logger.warning(
            "  kingfisher output for %s has %d invalid FASTQ(s) (%s); "
            "deleting all %d file(s) to force clean re-download",
            srr,
            len(bad),
            ", ".join(p.name for p in bad),
            len(produced),
        )
        for fq in produced:
            fq.unlink(missing_ok=True)
        return False
    return True


# ---------------------------------------------------------------------------
# Backend: sra-toolkit (optional)
# ---------------------------------------------------------------------------


def download_srr_sratoolkit(srr: str, output_dir: Path, dry_run: bool = False) -> bool:
    """Download an SRR via NCBI sra-toolkit (``prefetch`` + ``fasterq-dump``)."""
    if _skip_if_complete_else_clean(srr, output_dir):
        return True

    logger.info("  sra-toolkit %s → %s", srr, output_dir)
    if dry_run:
        return True

    try:
        subprocess.run(
            ["prefetch", srr, "-O", str(output_dir)],
            check=True,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        sra_cache = _find_sra_cache(output_dir, srr)
        fasterq_input = str(sra_cache) if sra_cache else srr
        if sra_cache is None:
            logger.warning(
                "  prefetch completed for %s but no local .sra under %s; "
                "trying fasterq-dump via accession",
                srr,
                output_dir,
            )

        subprocess.run(
            ["fasterq-dump", fasterq_input, "-O", str(output_dir), "--skip-technical", "--split-3"],
            check=True,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        for fq in iter_srr_files(output_dir, srr, suffix=".fastq"):
            subprocess.run(["gzip", str(fq)], check=True)
        if sra_cache is not None and sra_cache.exists():
            sra_cache.unlink()
    except subprocess.CalledProcessError as e:
        logger.error(
            "  sra-toolkit failed for %s: %s",
            srr,
            e.stderr[-300:] if e.stderr else str(e),
        )
        return False

    produced = iter_srr_files(output_dir, srr)
    if not produced:
        logger.warning("  sra-toolkit finished for %s but no FASTQ produced", srr)
        return False
    bad = [fq for fq in produced if not validate_fastq_gz(fq)]
    if bad:
        logger.warning(
            "  sra-toolkit output for %s has %d invalid FASTQ(s) (%s); deleting all %d file(s)",
            srr,
            len(bad),
            ", ".join(p.name for p in bad),
            len(produced),
        )
        for fq in produced:
            fq.unlink(missing_ok=True)
        return False
    return True


# ---------------------------------------------------------------------------
# Dispatcher (ENA-direct first, then optional kingfisher / sra-toolkit)
# ---------------------------------------------------------------------------


_KINGFISHER_GPL_WARNING = (
    "Using kingfisher backend (GPL-3.0). kingfisher is NOT bundled with "
    "selexprep's MIT install; it was invoked because it is present on PATH. "
    "If you redistribute selexprep outputs alongside kingfisher binaries you "
    "must respect GPL-3.0 terms."
)


def download_srr(
    srr: str,
    output_dir: Path,
    dry_run: bool = False,
    backend: DownloadBackend = "auto",
) -> bool:
    """Download one SRR via the selected backend.

    Parameters
    ----------
    backend
        - ``"auto"`` (default): ENA-direct → kingfisher (if installed) →
          sra-toolkit. MIT-compatible without external tools.
        - ``"ena"``: ENA-direct only. Fail fast if ENA can't serve. Use
          this in CI / paper benchmarks where GPL fallback would change
          the experiment.
        - ``"kingfisher"``: force kingfisher. **Emits a GPL-3.0 license
          notice via WARNING log.** Requires the binary on PATH.
        - ``"sra"``: force sra-toolkit (``prefetch`` + ``fasterq-dump``).

    Returns True on success, False after exhausting available backends.
    """
    if backend == "ena":
        logger.info("  ENA-direct (forced) for %s", srr)
        return download_srr_ena_direct(srr, output_dir, dry_run)

    if backend == "kingfisher":
        logger.warning(_KINGFISHER_GPL_WARNING)
        if not kingfisher_available():
            logger.error("  kingfisher backend requested but binary not on PATH")
            return False
        return download_srr_kingfisher(srr, output_dir, dry_run)

    if backend == "sra":
        if not sratoolkit_available():
            logger.error("  sra backend requested but prefetch+fasterq-dump not on PATH")
            return False
        return download_srr_sratoolkit(srr, output_dir, dry_run)

    # backend == "auto"
    logger.info("  trying ENA-direct for %s", srr)
    if download_srr_ena_direct(srr, output_dir, dry_run):
        return True
    logger.warning("  ENA-direct failed for %s", srr)

    if kingfisher_available():
        logger.warning(_KINGFISHER_GPL_WARNING)
        logger.info("  trying kingfisher fallback for %s", srr)
        if download_srr_kingfisher(srr, output_dir, dry_run):
            return True
        logger.warning("  kingfisher failed for %s", srr)

    if sratoolkit_available():
        logger.info("  trying sra-toolkit fallback for %s", srr)
        if download_srr_sratoolkit(srr, output_dir, dry_run):
            return True
        logger.warning("  sra-toolkit failed for %s", srr)

    logger.error("  all backends exhausted for %s", srr)
    return False


# ---------------------------------------------------------------------------
# External processed data (Zenodo / Figshare)
# ---------------------------------------------------------------------------


def zenodo_expected_md5(checksum: str | None) -> str | None:
    """Parse a Zenodo checksum string (``md5:hex...``) → bare hex, or None."""
    if not checksum:
        return None
    if checksum.startswith("md5:"):
        return checksum[4:]
    if len(checksum) == 32 and all(c in "0123456789abcdefABCDEF" for c in checksum):
        return checksum
    return None


def download_external_processed(
    bioproject_id: str,
    source: str,
    bp_raw_dir: Path,
    dry_run: bool = False,
) -> None:
    """Download Zenodo or Figshare processed-data files for a BioProject."""
    out_dir = bp_raw_dir / "processed_external"

    if source.startswith("zenodo:"):
        record_id = source.split(":")[1]
        url = f"https://zenodo.org/api/records/{record_id}"
        if dry_run:
            return
        try:
            meta = requests.get(url, timeout=30).json()
        except requests.RequestException as e:
            logger.warning("    Zenodo fetch failed for %s: %s", record_id, e)
            return
        if not meta:
            return
        for f in meta.get("files", []):
            file_url = f["links"]["self"]
            fname = f["key"]
            dest = out_dir / fname
            expected_md5 = zenodo_expected_md5(f.get("checksum"))
            if dest.exists() and expected_md5 and md5_file(dest) == expected_md5:
                logger.info("    [skip] %s already downloaded (md5 ok)", fname)
                continue
            logger.info("    downloading Zenodo file: %s", fname)
            stream_download(file_url, dest, expected_md5=expected_md5)

    elif source.startswith("figshare:"):
        article_id = source.split(":")[1]
        url = f"https://api.figshare.com/v2/articles/{article_id}/files"
        if dry_run:
            return
        try:
            files = requests.get(url, timeout=30).json()
        except requests.RequestException as e:
            logger.warning("    Figshare fetch failed for %s: %s", article_id, e)
            return
        for f in files:
            fname = f["name"]
            dest = out_dir / fname
            expected_md5 = f.get("computed_md5") or f.get("supplied_md5") or None
            if dest.exists() and expected_md5 and md5_file(dest) == expected_md5:
                logger.info("    [skip] %s already downloaded (md5 ok)", fname)
                continue
            logger.info("    downloading Figshare file: %s", fname)
            stream_download(f["download_url"], dest, expected_md5=expected_md5)


# ---------------------------------------------------------------------------
# Per-BioProject orchestrator
# ---------------------------------------------------------------------------


def _missing_srrs(raw_root: Path, bp_id: str, samples: list[dict]) -> list[str]:
    missing = []
    for sample in samples:
        srr = (sample.get("srr") or "").strip()
        if srr and not srr_present(raw_root, bp_id, srr):
            missing.append(srr)
    return missing


def download_bioproject(
    bp: dict,
    samples: list[dict],
    rounds_by_srr: dict[str, dict],
    raw_root: Path,
    processed_root: Path | None = None,
    coverage_info: dict | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Download every SRR for a BioProject.

    `raw_root` is the parent directory under which per-BP raw FASTQ
    subdirectories are written. `processed_root` is optional and used only
    for the resume-skip check (``summary.json`` presence). Post-processing
    (counting / enrichment) is NOT triggered here — that's a caller concern.
    """
    bp_id = bp["bioproject_id"]
    source = bp.get("source", "")

    status: dict[str, Any] = {
        "bioproject_id": bp_id,
        "n_srr": 0,
        "n_downloaded": 0,
        "n_failed": 0,
        "skipped": False,
    }

    if processed_root is not None and not dry_run:
        summary = processed_root / safe_dir_name(bp_id) / "summary.json"
        if summary.exists():
            missing = _missing_srrs(raw_root, bp_id, samples)
            if not missing:
                logger.info("[%s] already processed and all SRRs present — skipping", bp_id)
                status["skipped"] = True
                return status
            logger.warning(
                "[%s] summary.json exists but %d SRR(s) still missing; resuming download",
                bp_id,
                len(missing),
            )

    bp_raw_dir = raw_root / safe_dir_name(bp_id)
    bp_raw_dir.mkdir(parents=True, exist_ok=True)

    info_path = bp_raw_dir / "study_info.json"
    if not info_path.exists() and not dry_run:
        study_info = dict(bp)
        if coverage_info:
            study_info["round_coverage_status"] = coverage_info.get("round_coverage_status")
            study_info["round_coverage"] = coverage_info
        info_path.write_text(json.dumps(study_info, indent=2, ensure_ascii=False), encoding="utf-8")

    if source.startswith(("zenodo:", "figshare:")):
        download_external_processed(bp_id, source, bp_raw_dir, dry_run)
        return status

    if not samples:
        logger.warning("[%s] no samples provided", bp_id)
        return status

    status["n_srr"] = len(samples)

    for sample in samples:
        srr = sample.get("srr", "")
        if not srr:
            continue

        round_rec = rounds_by_srr.get(srr, {})
        rn_str = round_rec.get("round_number", "")
        rn = int(rn_str) if rn_str.strip().lstrip("-").isdigit() else None
        if needs_manual_review(round_rec):
            if rn is None:
                logger.warning("  [%s] round unknown — downloading to round_unknown/", srr)
            else:
                logger.warning(
                    "  [%s] round %02d assigned, but metadata needs manual review", srr, rn
                )
        target_hint = (sample.get("target_hint") or "").strip() or None

        out_dir = round_dir(raw_root, bp_id, rn, target_hint)
        ok = download_srr(srr, out_dir, dry_run)

        if ok:
            status["n_downloaded"] += 1
        else:
            status["n_failed"] += 1

    logger.info(
        "[%s] downloaded %d/%d SRRs (%d failed)",
        bp_id,
        status["n_downloaded"],
        status["n_srr"],
        status["n_failed"],
    )
    return status


# ---------------------------------------------------------------------------
# Top-level batch entry
# ---------------------------------------------------------------------------


def run_download(
    bioprojects: list[dict],
    samples_by_bp: dict[str, list[dict]],
    rounds_by_srr: dict[str, dict],
    raw_root: Path,
    processed_root: Path | None = None,
    coverage_filter: Callable[[dict], dict | None] | None = None,
    dry_run: bool = False,
) -> list[dict]:
    """Download every BioProject in `bioprojects`.

    `bioprojects` should already be filtered to those the caller wants to
    download (e.g. ``include == 'y'``). `samples_by_bp` and `rounds_by_srr`
    are the indexed lookups produced from a discovery run.

    `coverage_filter` is an optional callable returning the per-BP coverage
    dict (used in study_info.json); pass `None` to skip coverage annotation.
    """
    logger.info("Downloading %d BioProject(s)%s", len(bioprojects), " [DRY RUN]" if dry_run else "")

    statuses: list[dict] = []
    for bp in bioprojects:
        bp_id = bp["bioproject_id"]
        coverage = coverage_filter(bp) if coverage_filter else None
        logger.info(
            "=== %s: %s%s ===",
            bp_id,
            bp.get("protein_target", "?"),
            (f" [{coverage.get('round_coverage_status', 'unknown')}]" if coverage else ""),
        )
        bp_samples = samples_by_bp.get(bp_id, [])
        st = download_bioproject(
            bp,
            bp_samples,
            rounds_by_srr,
            raw_root=raw_root,
            processed_root=processed_root,
            coverage_info=coverage,
            dry_run=dry_run,
        )
        statuses.append(st)
        time.sleep(1)

    total_srr = sum(s["n_srr"] for s in statuses)
    total_ok = sum(s["n_downloaded"] for s in statuses)
    total_fail = sum(s["n_failed"] for s in statuses)
    skipped = sum(1 for s in statuses if s["skipped"])

    logger.info(
        "=== DOWNLOAD COMPLETE: %d BioProjects (%d skipped), %d/%d SRRs downloaded, %d failed ===",
        len(statuses),
        skipped,
        total_ok,
        total_srr,
        total_fail,
    )
    return statuses
