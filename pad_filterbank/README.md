# Pad SIGPROC filterbank snippets for folding #

This directory contains `python` code to pad a `SIGPROC` filterbank file to a certain length so that it can be folded with `DSPSR`. This is useful for filterbank snippets that are shorter than the period of a given source and/or when no timing ephemeris exists.

The code tries hard to reduce the mismatch between the data and the padding in overall normalisation and spectral behaviour (bandpass).

Optionally, the code can run `IQRM` RFI excision methods on the input data before padding. That is because it is often beneficial to filter out RFI before folding the data, or any other processing.

## Requirements ##

* [IQRM](https://github.com/v-morello/iqrm) (optional but recommended)
* Scipy
* [Your Unified Reader -- Your](https://github.com/thepetabyteproject/your)

## Usage ##

```bash
$ python3 pad_filterbank.py -h
usage: pad_filterbank.py [-h] -l seconds [--iqrm] filename

Pad the sigproc filterbank data.

positional arguments:
  filename              The name of the input filterbank file.

optional arguments:
  -h, --help            show this help message and exit
  -l seconds, --length seconds
                        The length in seconds to pad the filterbank data to. (default: None)
  --iqrm                Enable IQRM RFI excision prior to padding. (default: False)


python3 pad_filterbank.py 2020_09_27_01\:24\:08.fil -l 13.7 --iqrm
```
