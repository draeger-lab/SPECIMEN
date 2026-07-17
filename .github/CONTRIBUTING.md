# SPECIMEN CONTRIBUTING.md

## Welcome

Thank you for helping improve SPECIMEN. SPECIMEN is a collection of workflows for automated, standardized curation of genome-scale metabolic models. Contributions are most useful when they keep complete workflows reproducible, configurable, and documented for other modelling projects.

This guide focuses on the kinds of changes contributors usually make in SPECIMEN:

- Workflows and pipeline modules
- Configuration templates and default parameters
- Documentation and pipeline diagrams
- Example datasets and example input/output
- Interfaces to external tools
- Validation of complete pipeline runs

## Ways to contribute

You can contribute by:

- Reporting reproducible bugs in CMPB, HQTB, or future workflows
- Suggesting new workflow steps or refinements
- Improving configuration templates in `src/specimen/data/config/`
- Updating Sphinx documentation in `docs/source/`
- Adding or improving example notebooks in `HowTo/`
- Documenting required software, versions, and computational requirements
- Improving interfaces to external tools such as refineGEMs, CarveMe, DIAMOND, EntrezDirect, or ModelPolisher
- Adding example input/output that helps others validate a workflow change

## Reporting bugs

Please open a bug report with enough information for maintainers to reproduce the problem. Include:

- The workflow and command you ran, for example CMPB or HQTB through `specimen`
- The configuration file used, with private paths or accessions removed if needed
- Input data type and a minimal example, when possible
- SPECIMEN version or commit, Python version, and operating system
- External tool versions and installation method, especially Docker, pip editable install, or system packages
- The expected output and the actual output
- Logs, tracebacks, and generated reports
- Whether the same input worked with an earlier SPECIMEN version or configuration template

For workflow bugs, please say whether the failure happens before, during, or after a complete pipeline run. Partial-step failures are useful, but complete-run context helps identify configuration compatibility and downstream effects.

## Suggesting new features

Feature requests should describe the workflow problem first, then the proposed implementation. Helpful requests include:

- The modelling or curation task the feature supports
- The workflow step where it belongs
- Required input files and expected output files
- New external software, databases, web services, or credentials
- Expected runtime, memory, storage, or network requirements
- Compatibility expectations for existing configuration templates
- Documentation, diagram, and example-data updates needed to make the feature usable

## Development workflow

1. Create a feature branch from the active development branch.
2. Install SPECIMEN in editable mode with the dependencies needed for the workflow you are changing.
3. Make a focused change to one workflow, configuration area, documentation area, or external-tool interface at a time.
4. Update or add configuration examples for any new option.
5. Run the smallest relevant checks while developing, then validate at least one complete affected pipeline run before opening a pull request.
6. Update documentation, example notebooks, and diagrams when behavior visible to users changes.
7. Open a pull request with clear validation notes and any known limitations.

Avoid changing defaults or configuration keys without a compatibility plan. Existing workflows and example configurations should continue to run unless the pull request explicitly documents a breaking change.

## Expected directory structure

Use the existing project layout when adding or changing files:

- `src/specimen/cmpb/` contains CMPB workflow code.
- `src/specimen/hqtb/` contains HQTB workflow code.
- `src/specimen/hqtb/core/` contains HQTB pipeline modules.
- `src/specimen/data/config/` contains packaged configuration templates.
- `docs/source/` contains Sphinx documentation.
- `docs/source/images/` contains pipeline diagrams and documentation images.
- `HowTo/` contains example notebooks.
- `dev/` contains developer notes, experiments, and local helper material.

If a new workflow is added, keep its command-line interface, configuration templates, documentation pages, and diagrams discoverable from the same top-level documentation structure as the existing workflows.

## Adding new workflow steps

When adding a workflow step or pipeline module:

- Define its inputs, outputs, and side effects.
- Document how it is enabled, disabled, or configured.
- Keep file and directory naming deterministic.
- Validate that the step can run as part of a complete pipeline, not only in isolation.
- Include clear behavior for missing optional inputs, failed external commands, and empty intermediate results.
- Explain how the step interacts with existing CMPB, HQTB, or future workflow stages.
- Add or update example input/output when the step changes user-visible artifacts.

## Coding and documentation standards

- Follow the existing Python style and formatting used in the repository.
- Keep workflow code readable and explicit about data paths, generated files, and external commands.
- Use type hints and reStructuredText docstrings for public functions.
- Prefer configuration-driven behavior over hard-coded local paths.
- Keep default configuration files runnable for documented workflows.
- Document new configuration keys, allowed values, defaults, and compatibility notes.
- Update Sphinx pages and HowTo notebooks when the user-facing workflow changes.
- Update pipeline diagrams when step order, data flow, or major outputs change.

## Documenting required software

Changes that add or modify external dependencies should document:

- Software name and purpose
- Minimum tested version
- Installation source or reference documentation
- Required databases, indexes, credentials, or network access
- Expected command-line availability
- Runtime, memory, disk, and temporary-file requirements
- Docker-specific notes, if the dependency behaves differently in containers

Please include this information in the relevant workflow documentation and mention it in the pull request.

## Testing and validation

Validation should match the risk of the change. For workflow changes, maintainers need evidence that complete pipeline behavior still works.

Useful validation includes:

- Running the affected workflow from the command line with a documented configuration file
- Running a minimal example dataset through the full pipeline
- Comparing key output files, generated reports, and logs against expected results
- Checking that existing configuration templates still load and remain backward compatible
- Building or previewing affected Sphinx documentation
- Testing Docker behavior when the change affects installation, paths, or external tools

If a full run is too expensive, document the reason and provide the strongest smaller validation you could run.

## Pull request checklist

Before requesting review, please confirm:

- [ ] The change has a focused scope and a clear motivation.
- [ ] Existing workflows and configuration templates are not broken unexpectedly.
- [ ] New or changed configuration keys are documented.
- [ ] External dependencies, versions, and computational requirements are documented.
- [ ] Example input/output or notebooks are updated when user-facing behavior changes.
- [ ] Pipeline diagrams are updated when workflow structure changes.
- [ ] A complete affected pipeline run was tested, or a limitation is explained.
- [ ] Logs, generated reports, or other validation evidence are summarized in the pull request.
- [ ] The documentation builds or affected documentation files were reviewed.
- [ ] Breaking changes are called out explicitly.

## Review process

Maintainers will review contributions for:

- Reproducibility of the workflow or example
- Compatibility with existing configuration templates
- Robust handling of external tools and generated files
- Documentation completeness
- Clarity of validation evidence
- Impact on existing CMPB, HQTB, and future workflow behavior

Review may request smaller pull requests, more complete validation, additional documentation, or compatibility adjustments before merge.

## Getting help

If you are unsure where a change belongs, open an issue before starting a large implementation. For questions, include the workflow, configuration file, expected inputs and outputs, external tools involved, and what you have already tried.

## Citation

If you use SPECIMEN in research, please cite the project using the citation information in the repository and the Zenodo DOI badge shown in the README:

https://doi.org/10.5281/zenodo.12723500
