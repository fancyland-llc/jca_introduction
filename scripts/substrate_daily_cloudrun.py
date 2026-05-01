"""Cloud Run function entry point for the daily substrate pipeline.

Deployed as a SEPARATE Cloud Run service (not merged into prompt-studio-backend).
Triggered daily by Cloud Scheduler via HTTP POST.

Architecture:
  Cloud Scheduler (cron: 02:00 UTC daily)
    -> HTTP POST -> this Cloud Run service
    -> runs substrate_daily.run_pipeline()
    -> writes daily_status doc + fiedler_history entry to Firestore
    -> customer dashboard reads from Firestore

Deployment (via gcloud, NOT merged into prompt-studio-backend):
  gcloud run deploy substrate-daily-pipeline \\
    --source . \\
    --function substrate_daily_pipeline \\
    --base-image python313 \\
    --region us-central1 \\
    --memory 2Gi \\
    --cpu 2 \\
    --timeout 1800 \\
    --no-allow-unauthenticated \\
    --service-account substrate-runner@brainstormingorganization.iam.gserviceaccount.com \\
    --project brainstormingorganization

Or via the Cloud Console "Write a function" path: select Python 3.13,
function name `substrate_daily_pipeline`, paste this file's contents +
substrate_daily.py + supporting modules, set memory 2Gi / timeout 30 min.

Cloud Scheduler trigger setup:
  gcloud scheduler jobs create http substrate-daily-trigger \\
    --schedule "0 2 * * *" \\
    --time-zone "UTC" \\
    --uri "https://substrate-daily-pipeline-XXX.run.app" \\
    --http-method POST \\
    --oidc-service-account-email substrate-scheduler@brainstormingorganization.iam.gserviceaccount.com
"""

from __future__ import annotations

import json
import logging
import traceback
from datetime import datetime, timezone

import functions_framework

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("substrate_daily_cloudrun")


@functions_framework.http
def substrate_daily_pipeline(request):
    """HTTP entry point for the daily substrate pipeline.

    Idempotent: safe to call multiple times per day; latest run wins in
    Firestore. Triggered daily by Cloud Scheduler via OIDC-authenticated
    POST. Manual trigger via gcloud:

      curl -X POST -H "Authorization: Bearer $(gcloud auth print-identity-token)" \\
        https://substrate-daily-pipeline-XXX.run.app
    """
    started_at = datetime.now(timezone.utc)
    log.info("substrate_daily_pipeline trigger received at %s", started_at.isoformat())

    # Optional request-body overrides (for manual / debug triggers)
    overrides: dict = {}
    try:
        if request.method == "POST" and request.is_json:
            overrides = request.get_json(silent=True) or {}
    except Exception:
        overrides = {}

    baseline_only = bool(overrides.get("baseline_only", False))
    write_firestore = bool(overrides.get("write_firestore", True))
    k_nn = int(overrides.get("k_nn", 80))

    try:
        # Late import: keeps cold-start cost off the trigger; the heavyweight
        # numpy/scipy/firebase_admin imports happen inside substrate_daily.
        from substrate_daily import run_pipeline
        status = run_pipeline(
            write_firestore=write_firestore,
            baseline_only=baseline_only,
            k_nn=k_nn,
            verbose=True,
        )
    except Exception as e:
        log.exception("substrate pipeline failed")
        return (
            json.dumps({
                "status": "error",
                "error": str(e),
                "traceback": traceback.format_exc(),
                "started_at": started_at.isoformat(),
            }),
            500,
            {"Content-Type": "application/json"},
        )

    finished_at = datetime.now(timezone.utc)
    duration_sec = (finished_at - started_at).total_seconds()
    log.info("substrate_daily_pipeline complete in %.1fs", duration_sec)

    # Compact response: full payload is in Firestore; HTTP response is for
    # Cloud Scheduler health checks + manual trigger debugging.
    response = {
        "status": "ok",
        "started_at": started_at.isoformat(),
        "finished_at": finished_at.isoformat(),
        "duration_sec": duration_sec,
        "lambda_2_baseline": status.get("fiedler_baseline", {}).get("lambda_2"),
        "n_nodes": status.get("fiedler_baseline", {}).get("n_nodes"),
        "n_edges": status.get("fiedler_baseline", {}).get("n_edges"),
        "psigma_leo": status.get("catalog", {}).get("psigma_leo"),
        "current_kp": status.get("noaa", {}).get("current_kp"),
        "stress_test_count": len(status.get("fiedler_stress_tests") or {}),
        "any_ballistic": any(
            v.get("ballistic_transition", False)
            for v in (status.get("fiedler_stress_tests") or {}).values()
        ),
        "firestore_doc": "substrate_daily/current",
    }
    return (json.dumps(response), 200, {"Content-Type": "application/json"})
