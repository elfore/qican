#!/usr/bin/env python3
"""Compatibility entry point for the controlled YY1 assessment workflow.

Use ``--apply-feishu`` only after the local 15-table CSV validation succeeds.
"""

from run_assessment_workflow_v2 import main


if __name__ == "__main__":
    raise SystemExit(main())
