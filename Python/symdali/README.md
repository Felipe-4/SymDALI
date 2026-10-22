# symdali

symdali is a Python extension written in C that provides functions for the [SymDALI](https://github.com/Felipe-4/SymDALI) framework.

## Requirements

- Python 3.11 or later
- NumPy 1.26 or later, but below 2.0
- A compatible C compiler

Installing from source also requires the Python development headers. Depending on your operating system, you may need to install additional build tools.

## Installation

Install symdali directly from PyPI:

```bash
python -m pip install symdali
```

symdali is distributed as source code. During installation, `pip` compiles the C extension on your machine.


## Usage

See [GitHub](https://github.com/Felipe-4/SymDALI) and links therein for the documentation.


## Platform support

Source compilation has been tested on:

- macOS ARM64 with Python 3.11
- Linux x86_64 with Python 3.12

Other platforms may work if the required Python development headers, NumPy version, and C compiler are available. Windows has not been tested.

## License

See the [LICENSE](LICENSE) file for the terms under which symdali is distributed.
