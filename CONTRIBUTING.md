# Contributing to MOPEDS

Thank you for your interest in contributing to MOPEDS.

Contributions of all sizes are welcome, including:

- Bug reports
- Feature requests
- Documentation improvements
- Examples and case studies
- Tests
- Code contributions

## Development Setup

Clone the repository:

```bash
git clone https://git.tu-berlin.de/dbta/optimization/mopeds.git
cd mopeds
```

Install development dependencies:

```bash
uv sync --all-groups
```

Run the test suite:

```bash
uv run pytest
```

Run tests with coverage or generate report:

```bash
uv run pytest --cov=mopeds --cov-report=term-missing
uv run pytest --cov=mopeds --cov-report=html
```


## Contribution Workflow

### GitLab

The primary repository of MOPEDS is hosted on TU Berlin GitLab.

Users with access to GitLab may create a branch and submit a merge request directly.

### GitHub

A public GitHub mirror is available for community contributions.

If you do not have a GitLab account, please:

1. Fork the GitHub repository.
2. Create a feature branch in your fork.
3. Submit a pull request.

Pull requests submitted through GitHub are reviewed and synchronized with the primary GitLab repository.

Please create a separate branch for your work and submit a merge request (GitLab) or pull request (GitHub).

Before submitting, ensure that:

- All tests pass.
- New functionality is documented.
- Public APIs include appropriate docstrings.
- Code follows the existing project style.

## Requirements for New Features

Every new feature should include:

- At least one automated test.
- At least one runnable example demonstrating the feature.
- Documentation updates if the public API changes.

Contributions that add functionality without tests or examples are unlikely to be accepted.

## Bug Fixes

Bug fixes should include a regression test whenever possible.

## Questions

If you are unsure how to implement a feature or fix, feel free to open an issue before starting work.

## Commit Message Style

MOPEDS uses [Conventional Commits](https://www.conventionalcommits.org/) to generate the changelog automatically.

Please write commit messages in the following format:

```text
type: short description
```

Common types are:

```text
feat: add a new feature
fix: fix a bug
docs: update documentation
test: add or update tests
refactor: change code without changing behavior
ci: update CI configuration
```

Breaking changes should be marked explicitly:

```text
feat!: require CasADi 3.7
```

These commit messages are used to generate `CHANGELOG.md`, so please keep them clear and user-facing.

To generate `CHANGELOG.md` run:

```bash
uv run cz changelog
```
