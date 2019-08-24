# Numerov

A basic implementation of the Numerov method to solve the Schrodinger equation for a 1/r potential

## Install

### download source

```
git clone "https://github.com/ad3ller/numerov"
cd ./numerov
```

### cython build (optional)

The Numerov method can be sped up significantly (x10) with [cython](https://cython.org/). (Requires a C compiler, e.g., gcc or Visual Studio SDK).

```
python setup.py build_ext -i
```

### install

```
python setup.py install
```

### run tests (optional)

Requires cython build, sympy and pytest.

```
pytest
```

## Docs

```
>>> import numerov
>>> numerov.radial_integral(12, 5, 15, 4, step=0.0001)
4.573187231242028
```

See notebooks.
