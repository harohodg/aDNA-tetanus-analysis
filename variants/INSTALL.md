To use these scripts, you will need one of the following virutalization software options:
- [docker](https://docs.docker.com/), or
- [singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)

When you have one of these installed, use `./scripts/setup.sh [docker/singularity]` to
initialize your work environment. This will pull the appropriate images.

The full analysis is wrapped in a single script, `variants.sh`, which can be run via either Docker or Singuarity.
