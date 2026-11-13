Contributing
============

Contributions to ``SPECIMEN`` are welcome. Because ``SPECIMEN`` is built around
complete modelling workflows, contributions usually focus on reproducibility and
workflow behaviour rather than only on Python APIs.

Typical contribution areas include:

* workflow definitions and pipeline modules,
* configuration templates and default parameters,
* documentation and pipeline diagrams,
* example datasets and example input/output,
* interfaces to external tools.

Community Files
---------------

The main contribution guidance is maintained in the GitHub community files:

* `Contributing guide <https://github.com/draeger-lab/SPECIMEN/blob/dev/.github/CONTRIBUTING.md>`__
* `Bug report template <https://github.com/draeger-lab/SPECIMEN/blob/dev/.github/ISSUE_TEMPLATE/bug_report.yml>`__
* `Feature request template <https://github.com/draeger-lab/SPECIMEN/blob/dev/.github/ISSUE_TEMPLATE/feature_request.yml>`__
* `Workflow suggestion template <https://github.com/draeger-lab/SPECIMEN/blob/dev/.github/ISSUE_TEMPLATE/workflow_suggestion.yml>`__
* `Pull request template <https://github.com/draeger-lab/SPECIMEN/blob/dev/.github/pull_request_template.md>`__
* `Code of Conduct <https://github.com/draeger-lab/SPECIMEN/blob/dev/.github/CODE_OF_CONDUCT.md>`__

Before opening a pull request, please check the contribution guide and use the
provided templates. These files help contributors and maintainers communicate
the information needed to review workflow changes consistently.

What to Document
----------------

For changes to workflows, configuration templates, examples, or external-tool
interfaces, please document:

* the affected workflow and pipeline step,
* required input files and generated output files,
* new or changed configuration keys,
* required software, databases, versions, and installation notes,
* expected runtime, memory, disk, and network requirements,
* example input/output, when possible,
* updates to Sphinx documentation, HowTo notebooks, or pipeline diagrams.

Validation Expectations
-----------------------

Workflow changes should be validated with a complete affected pipeline run when
possible. If a complete run is too expensive or requires unavailable data,
please describe the smaller validation that was performed and explain the
limitation in the pull request.

Useful validation includes checking that existing configuration templates still
load, expected output files are generated, external tools are called correctly,
and existing workflows are not broken unexpectedly.
