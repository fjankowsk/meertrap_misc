#
#   2021 Fabian Jankowski
#   Pad sigproc filterbank data to a given length.
#

import argparse
import os.path

import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import shape
from numpy.lib.arraypad import pad
from scipy import interpolate
from iqrm import iqrm_mask
import your
from your.formats.filwriter import make_sigproc_object


def parse_args():
    """
    Parse the commandline arguments.

    Returns
    -------
    args: populated namespace
        The commandline arguments.
    """

    parser = argparse.ArgumentParser(
        description='Pad the sigproc filterbank data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        'filename',
        type=str,
        help='The name of the input filterbank file.'
    )

    parser.add_argument(
        '-l',
        '--length',
        dest='length',
        metavar=('seconds'),
        type=float,
        required=True,
        help='The length in seconds to pad the filterbank data to.'
    )

    parser.add_argument(
        '--iqrm',
        action='store_true',
        dest='iqrm',
        default=False,
        help='Enable IQRM RFI excision prior to padding.'
    )

    args = parser.parse_args()

    return args


def plot_bandpass_data(bandpass, spline, padding_data):
    """
    Plot the bandpass data.

    Parameters
    ----------
    bandpass: ~np.array
        The bandpass data.
    spline: ~scipy.interpolate.UnivariateSpline
        The best spline fit to the data.
    """

    samples = np.arange(len(bandpass))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(
        bandpass,
        color='C0',
        lw=1,
        label='data',
        zorder=3
    )

    ax.plot(
        samples,
        spline(samples),
        color='C1',
        lw=1,
        ls='dashed',
        label='spline fit',
        zorder=4
    )

    ax.plot(
        padding_data,
        color='C2',
        lw=1,
        ls='dotted',
        label='padding data',
        zorder=4
    )

    ax.set_xlabel('Channel number')
    ax.set_ylabel('Amplitude')
    ax.set_title('Median bandpass')
    ax.legend(loc='best', frameon=False)

    fig.tight_layout()


def pad_data(filename, length, iqrm):
    """
    Pad the data as necessary to reach a certain length.

    Parameters
    ----------
    filename: str
        The name of the input sigproc filterbank file.
    length: float
        The length in seconds to pad the data to.
    iqrm: bool
        Whether to run IQRM RFI excision.
    """

    yobj = your.Your(filename)
    print(yobj.your_header)

    data = yobj.get_data(
        nstart=0,
        nsamp=yobj.your_header.nspectra
    )

    # run iqrm rfi excision
    if iqrm:
        mean = np.mean(data, axis=None)
        spectral_std = np.std(data, axis=0)
        mask, _ = iqrm_mask(spectral_std, radius=2)
        print('IQRM channel mask: {}'.format(np.where(mask)[0]))

        for isamp in range(data.shape[0]):
            data[isamp, mask[0]] = mean

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(
        np.mean(data, axis=1),
        lw=0.5,
        color='black'
    )

    fig.tight_layout()

    mean = np.mean(data, axis=None)
    median = np.median(data, axis=None)
    quantiles = np.quantile(data, q=[0.25, 0.75], axis=None)
    std = 0.7413 * np.abs(quantiles[1] - quantiles[0])
    print('Mean, median and robust std of input filterbank: {0:.4f}, {1:.4f}, {2:.4f}'.format(mean, median, std))

    # fit median bandpass
    bandpass = np.median(data, axis=0)
    quantiles = np.quantile(data, q=[0.25, 0.75], axis=0)
    bandpass_std = 0.7413 * np.abs(quantiles[1] - quantiles[0])
    samples = np.arange(len(bandpass))

    mask = (bandpass >= median - 2 * std) & \
           (bandpass <= median + 2 * std)

    spline = interpolate.UnivariateSpline(
        x=samples[mask],
        y=bandpass[mask],
        k=2,
        s=20 * len(bandpass)
    )

    nspectra_padding = int(np.ceil(length / yobj.your_header.tsamp - yobj.your_header.nspectra))
    print('Number of spectra to append: {0}'.format(nspectra_padding))

    padding_data = np.zeros(shape=(nspectra_padding, yobj.your_header.nchans))
    noise = np.zeros(shape=(nspectra_padding, yobj.your_header.nchans))

    rng = np.random.default_rng()
    for ichan in range(yobj.your_header.nchans):
        noise[:, ichan] = rng.normal(
            loc=0,
            scale=bandpass_std[ichan],
            size=nspectra_padding
        )

    # apply the bandpass model to the padding data
    for isamp in range(nspectra_padding):
        for ichan in range(yobj.your_header.nchans):
            padding_data[isamp, ichan] = bandpass[ichan] + noise[isamp, ichan]

    plot_bandpass_data(bandpass, spline, np.median(padding_data, axis=0))

    padding_data = np.round(padding_data)
    padding_data = padding_data.astype(yobj.your_header.dtype)

    print('Shape of padding data: {0}'.format(padding_data.shape))

    # centre the data in amplitude and match the means
    shift = int(np.round(128 - mean))
    data = data + shift

    shift = int(np.round(128 - np.mean(padding_data, axis=None)))
    padding_data = padding_data + shift

    outname = '{0}_padded.fil'.format(
        os.path.splitext(filename)[0]
    )

    sigproc = make_sigproc_object(
        rawdatafile=outname,
        source_name=yobj.source_name,
        nchans=yobj.nchans,
        foff=yobj.foff,
        fch1=yobj.fch1,
        tsamp=yobj.your_header.tsamp,
        tstart=yobj.your_header.tstart,
        src_raj=yobj.your_header.ra_deg,
        src_dej=yobj.your_header.dec_deg,
        machine_id=0,
        nbeams=yobj.nbeams,
        ibeam=yobj.ibeam,
        nbits=yobj.your_header.nbits,
        nifs=1,
        barycentric=0,
        pulsarcentric=0,
        telescope_id=0,
        data_type=0,
        az_start=-1,
        za_start=-1,
    )

    sigproc.write_header(outname)

    sigproc.append_spectra(
        data,
        outname
    )

    sigproc.append_spectra(
        padding_data,
        outname
    )

    # check that all went fine
    padded_yobj = your.Your(outname)
    print(padded_yobj.your_header)

    assert (padded_yobj.your_header.nspectra == nspectra_padding + yobj.your_header.nspectra)

    for field in ['bw', 'dtype', 'nchans', 'tsamp', 'tstart']:
        assert (getattr(padded_yobj.your_header, field) == getattr(yobj.your_header, field))


#
# MAIN
#

def main():
    args = parse_args()

    pad_data(args.filename, args.length, args.iqrm)

    plt.show()

    print('All done.')


if __name__ == '__main__':
    main()
