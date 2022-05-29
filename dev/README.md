### Development documentation

#### Packaging
Linux wheels are built using the maturin manylinux docker container with:

```bash
docker run -it --rm -v $(pwd):/io --entrypoint bash ghcr.io/pyo3/maturin:main
maturin build --release --strip --interpreter python3.7 python3.8 python3.9 python3.10
```

while macOS wheels are generated with:
```bash
./dev/build_mac_wheels.sh
```

Upload to TestPyPI first with
```bash
python3 -m twine upload --repository testpypi target/wheels/*
```

