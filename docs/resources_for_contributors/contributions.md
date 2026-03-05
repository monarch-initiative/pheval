# Contributing to PhEval

Thank you for your interest in contributing to **PhEval**.  
Contributions are welcome across code, documentation, testing, and ecosystem plugins.

This page outlines **how to contribute**, **what kinds of contributions are encouraged**, and **how to get started**.

---

## Ways to contribute

You can contribute to PhEval in several ways:

- Reporting bugs or unexpected behaviour
- Suggesting enhancements or new features
- Improving documentation or examples
- Adding tests or improving coverage
- Developing new plugins or runners
- Maintaining or extending existing plugins

Not all contributions need to involve code - documentation and feedback are equally valuable.

---

## Before you start

Before opening an issue or pull request, please:

1. Check existing issues to avoid duplication
2. Ensure you are using a supported version of PhEval
3. Read the relevant documentation section (especially for plugins and runners)

If you are unsure whether something is a bug or a usage issue, opening a discussion or issue is encouraged.

---

## Reporting issues

Bugs, questions, and feature requests should be reported via GitHub issues.

When reporting an issue, please include:

- PhEval version
- Plugin name and version (if applicable)
- Command(s) executed
- Relevant configuration files (sanitised if needed)
- Error messages or stack traces
- Expected vs observed behaviour

Clear, minimal examples help issues get resolved faster.

---

## Contributing code

### General guidelines

- Keep changes focused and scoped
- Follow existing code style and patterns
- Add tests where appropriate
- Update documentation if behaviour changes

For larger changes, consider opening an issue first to discuss design and scope.

---

### Development setup

PhEval uses modern Python tooling. A typical development setup involves:

- Python 3.10+
- `uv` for dependency management
- `coverage` for testing
- `ruff` for linting and formatting

After cloning the repository:

```bash
uv sync
uv run coverage run -p -m pytest --durations=20 tests
```

Ensure all tests pass before submitting a pull request.

---

## Pull requests

When submitting a pull request:

- Clearly describe what the change does
- Reference any related issues
- Explain breaking changes explicitly
- Ensure CI passes

Pull requests should be small enough to review easily whenever possible.

---

## Plugin contributions

If you are contributing a new plugin:

- Use the PhEval runner template
- Follow the standard runner interface
- Ensure standardised result schemas are respected
- Verify result file stem matching with phenopackets
- Document tool-specific configuration clearly in the plugin README

See:

- [`Developing a PhEval plugin`](developing_a_pheval_plugin.md) for detailed guidance

---

## Documentation contributions

Documentation is built using **MkDocs Material**.

You can contribute by:

- Fixing typos or clarifying wording
- Adding examples
- Improving structure or navigation

Documentation-only pull requests are welcome.

---

## Getting help

If you are unsure where to start:

- Open an issue describing what you would like to contribute
- Ask questions in an existing issue or discussion
- Start with documentation improvements to familiarise yourself with the project

We appreciate all contributions, big or small, thank you for helping improve PhEval.
