# Contributing

Please feel free to contribute to `nisarhdf` by submitting issues or pull requests on GitHub.


## Building a new package version

1. Bump version tag and add additional dependencies if necessary in `pyproject.toml`

1. Push tag to GitHub that matches new version, for example:
```
git tag v0.4.1
git push --tags
```

1. Build the conda-package with pixi:
```
pixi build
```

1. Upload new package to prefix "channel" to install via `conda` tooling, for example:
```
pixi upload https://prefix.dev/api/v1/upload/fastice nisarhdf-0.4.1-pyh4616a5c_0.conda
```

NOTE: you must be an 'Owner' of the https://prefix.dev/channels/fastice channel to upload packages.

