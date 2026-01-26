# Release Guide for sp-pharkin

This document outlines the process for tagging and releasing new versions of sp-pharkin.

## Pre-Release Checklist

- [ ] All tests passing: `uv run pytest tests/ --cov=sp_pharkin -q`
- [ ] Code linted: `uv run ruff check . --fix`
- [ ] All code formatted: `uv run ruff format .`
- [ ] No Pylance errors: Review `.github/copilot-instructions.md` known issues
- [ ] CHANGELOG updated with version and features
- [ ] README updated if needed (API changes, new functions)
- [ ] All changes ready to include in release (no stashing)

## Version Numbering

Use [Semantic Versioning](https://semver.org/) (MAJOR.MINOR.PATCH):

- **MAJOR**: Breaking API changes (e.g., function signature changes)
- **MINOR**: New features, backward compatible (e.g., new functions)
- **PATCH**: Bug fixes, documentation updates, backward compatible

Current version: Check `pyproject.toml` under `[project] version = "X.Y.Z"`

## Release Steps

### 0. Review and Stage All Changes

First, review all uncommitted changes that will be included in the release:

```bash
git status
git diff
```

Stage all changes:

```bash
git add .
```

Verify staged changes:

```bash
git diff --cached
```

### 1. Update Version

Edit `pyproject.toml`:
```toml
[project]
version = "0.0.5"  # Increment version
```

### 2. Create Commit with Change Description

Build a descriptive commit message based on the staged changes. Include what was added/fixed/improved:

```bash
git commit -m "Release: v0.0.5

- New features and improvements
- Bug fixes
- Documentation updates"
```

Alternatively, let git open your editor for a more detailed message:

```bash
git commit
```

Add a detailed description of all changes (the staged diff should inform this message).

### 3. Create Git Tag

```bash
git tag -a v0.0.5 -m "Version 0.0.5"
git push origin master
git push origin v0.0.5
```

### 4. Create GitHub Release

On GitHub:
1. Go to Releases → Draft a new release
2. Select the tag you just created (v0.0.5)
3. Auto-generate release notes or manually write highlights
4. Publish release

## References

- [Semantic Versioning](https://semver.org/)
- [Git Tagging](https://git-scm.com/book/en/v2/Git-Basics-Tagging)

