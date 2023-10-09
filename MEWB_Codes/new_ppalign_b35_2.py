#!/usr/bin/env python

###########
# ppalign #
###########

# ppalign is a command-line program used to align homogeneous data (i.e. from
#    the same receiver, with the same center frequency, bandwidth, and number
#    of channels).  This is useful for making averaged portraits to either pass
#    to ppgauss.py with -M to make a Gaussian model, or to smooth and use as a
#    model with pptoas.py.

# Written by Timothy T. Pennucci (TTP; tim.pennucci@nanograv.org).

# Need option for constant Gaussian initial guess.

from __future__ import division
from __future__ import print_function

import os
import shlex
import subprocess as sub
from builtins import map
from builtins import range

from past.utils import old_div
from pptoas import *


def psradd_archives(metafile, outfile, palign=False):
    """
    Add together archives using psradd.

    This function will call psradd with an option to pass -P and can be used to
    make an initial guess for align_archives.

    metafile is a file containing PSRFITS archive names to be averaged.
    outfile is the name of the output archive.
    palign=True passes -P to psradd, which phase-aligns the archives, intead of
        using the ephemeris (maybe?).
    """
    psradd_cmd = "psradd "
    if palign:
        psradd_cmd += "-P "
    psradd_cmd += "-T -o %s -M %s" % (outfile, metafile)
    psradd_call = sub.Popen(shlex.split(psradd_cmd))
    psradd_call.wait()


def psrsmooth_archive(archive, options="-W"):
    """
    Smooth an archive using psrsmooth.

    This function will call psrsmooth with options to smooth an output archive
    from align_archives.

    archive is the PSRFITS archive to be smoothed.
    options are the options passed to psrsmooth.
    """
    psrsmooth_cmd = "psrsmooth " + options + " %s" % archive
    psrsmooth_call = sub.Popen(shlex.split(psrsmooth_cmd))
    psrsmooth_call.wait()


def align_archives(metafile, initial_guess, phase_offset, nu_mean_band, fit_dm=True, tscrunch=False,
        pscrunch=True, SNR_cutoff=0.0, outfile=None, norm=None, rot_phase=0.0,
        place=None, niter=1, quiet=False):
    """
    Iteratively align and average archives.

    Each archive is fitted for a phase, a DM, and channel amplitudes against
    initial_guess.  The average is weighted by the fitted channel amplitudes
    and channel noise level.  The average becomes the new initial alignment
    template for additional iterations.  The output archive will have a 0 DM
    value and dmc=0.

    metafile is a file containing PSRFITS archive names to be averaged, or it
        can be a python list of archive names.
    initial_guess is the PSRFITS archive providing the initial alignment guess.
    fit_dm=False will only align the subintegrations with a fit for phase.
    tscrunch=True will pre-average the subintegrations; recommended unless
        there is a reason to keep the invidual subints for looping over.
    pscrunch=False will return the average Stokes portraits.  Alignment and
        weighting is performed only via the total intensity portrait.
    SNR_cutoff is used to filter out low S/N archives from the average.
    outfile is the name of the output archive; defaults to
        <metafile>.algnd.fits.
    norm is the normalization method (None, 'mean', 'max', 'prof', 'rms', or
        'abs') applied to the final data.
    rot_phase is an overall rotation to be applied to the final output archive.
    place is a phase value at which to roughly place the peak pulse; it
        overrides rot_phase.
    niter is the number of iterations to complete.  1-5 seems to work ok.
    quiet=True suppresses output.

    """
    if type(metafile) == str:  # this was the old default
        datafiles = \
            [datafile[:-1] for datafile in open(metafile, "r").readlines()]
        if outfile is None:
            outfile = metafile + ".algnd.fits"
    else:  # assume metafile is list or array already; to be updated
        datafiles = metafile
    vap_cmd = "vap -c nchan,nbin %s" % initial_guess
    nchan, nbin = list(map(int, sub.Popen(shlex.split(vap_cmd), stdout=sub.PIPE
                                          ).stdout.readlines()[1].split()[-2:]))
    if pscrunch:
        state = 'Intensity'
        npol = 1
    else:
        state = 'Stokes'
        npol = 4
    try:
        model_data = load_data(initial_guess, state=state, dedisperse=False,
                               dededisperse=False, tscrunch=True, pscrunch=pscrunch,
                               fscrunch=False, rm_baseline=True, flux_prof=False,
                               refresh_arch=True, return_arch=True, quiet=quiet)
        # print(model_data)
        # plt.plot(figsize=(20,20))
        # plt.imshow(model_data,extent=(0,1,0,32), aspect='auto')
        # plt.title("portx_34")
        # plt.xlabel('nbin')
        # plt.ylabel('nchan')
        # plt.show()
        
    except IndexError:
        print("%s: has npol = 1; need npol == 4." % initial_guess)
        sys.exit()
    model_port_temp = (model_data.masks * model_data.subints)[0, 0]
    P_band = model_data.Ps[0]
    freqsx_band = model_data.freqs[0]
    DM = model_data.DM
    
    #print(np.shape(freqsx_band))
    #plt.plot(figsize=(20,20))
    #plt.imshow(model_port_temp,extent=(0,1,0,1), aspect='auto')
    #plt.title("portx_34")
    #plt.xlabel('nbin')
    #plt.ylabel('nchan')
    #plt.show()
    
    model_port = rotate_data(model_port_temp, phase_offset, DM/model_data.doppler_factors, P_band, freqsx_band,
                           nu_mean_band)  # why not nu_fit?
    #print(np.shape(model_port_temp))
    #plt.plot(figsize=(20,20))
    #plt.imshow(model_port,extent=(0,1,0,1), aspect='auto')
    #plt.title("portx_34")
    #plt.xlabel('nbin')
    #plt.ylabel('nchan')
    #plt.show()
    
    
    
    skip_these = []
    count = 1
    while (niter):
        print("Doing iteration %d..." % count)
        load_quiet = quiet
        aligned_port = np.zeros((npol, nchan, nbin))
        total_weights = np.zeros((nchan, nbin))
        if count == 2:
            for skipfile in skip_these:
                skipped = datafiles.pop(datafiles.index(skipfile))
        for ifile in range(len(datafiles)):
            try:
                data = load_data(datafiles[ifile], state=state,
                        dedisperse=False, tscrunch=tscrunch, pscrunch=pscrunch,
                        fscrunch=False, rm_baseline=rm_baseline,
                        flux_prof=False, refresh_arch=False, return_arch=False,
                        quiet=load_quiet)

            except RuntimeError:
                if not quiet:
                    print("%s: cannot load_data().  Skipping it." % \
                          datafiles[ifile])
                skip_these.append(datafiles[ifile])
                continue
            except IndexError:
                if not quiet:
                    print("%s: has npol = 1.  Skipping it." % \
                          datafiles[ifile])
                skip_these.append(datafiles[ifile])
                continue
            if data.nbin != model_data.nbin:
                if not quiet:
                    print("%s: %d != %d phase bins.  Skipping it." % \
                          (datafiles[ifile], data.nbin, model_data.nbin))
                skip_these.append(datafiles[ifile])
                continue
            if data.prof_SNR < SNR_cutoff:
                if not quiet:
                    print("%s: %d < %d S/N cutoff.  Skipping it." % \
                          (datafiles[ifile], data.prof_SNR, SNR_cutoff))
                skip_these.append(datafiles[ifile])
                continue
            try:
                freq_diffs = data.freqs - model_data.freqs
                if freq_diffs.min() == freq_diffs.max() == 0.0:
                    same_freqs = True
                else:
                    same_freqs = False
            except:
                same_freqs = False
            DM_guess = data.DM * np.logical_not(data.dmc)  # = 0.0 if de-disp'd
            for isub in data.ok_isubs:
                if same_freqs:
                    ichans = np.intersect1d(data.ok_ichans[isub],
                                            model_data.ok_ichans[0])
                    model_ichans = ichans
                else:  # Find 'close' frequency indices
                    ichans = data.ok_ichans[isub]
                    model_ichans = []
                    for ichan in ichans:
                        data_freq = data.freqs[isub, ichan]
                        imin = np.argmin(abs(model_data.freqs[0] - data_freq))
                        model_ichans.append(imin)
                    model_ichans = np.array(model_ichans)
                port = data.subints[isub, 0, ichans]
                freqs = data.freqs[isub, ichans]  # Use data freqs
                # freqs = model_data.freqs[isub,model_ichans]  #Use model freqs
                model = model_port[model_ichans]
                P = data.Ps[isub]
                SNRs = data.SNRs[isub, 0, ichans]
                errs = data.noise_stds[isub, 0, ichans]
                nu_fit = guess_fit_freq(freqs, SNRs)
                rot_port = rotate_data(port, 0.0, DM_guess, P, freqs,
                                       nu_fit)
                phase_guess = fit_phase_shift(np.average(rot_port, axis=0,
                                                         weights=data.weights[isub, ichans]), model.mean(axis=0),
                                              Ns=nbin).phase
                if len(freqs) > 1:
                    fit_phase = 1
                    fit_dm = int(bool(fit_dm))
                    fit_flags = [fit_phase, fit_dm, 0, 0, 0]
                    results = fit_portrait_full(port, model,
                            [phase_guess, DM_guess, 0.0, 0.0, 0.0], P, freqs,
                            [nu_fit, nu_fit, nu_fit], [None, None, None], errs,
                            fit_flags, log10_tau=False, quiet=quiet)
                    results.phase = results.phi
                    results.nu_ref = results.nu_DM
                else:  #1-channel hack
                    results = fit_phase_shift(port[0], model[0], errs[0],
                                              Ns=nbin)
                    results.DM = data.DM
                    results.nu_ref = freqs[0]
                    results.scales = np.array([results.scale])
                weights = np.outer(old_div(results.scales, errs ** 2), np.ones(nbin))
                for ipol in range(npol):
                    aligned_port[ipol, model_ichans] += weights * \
                            rotate_data(data.subints[isub,ipol,ichans],
                                    results.phase, results.DM, P, freqs,
                                    results.nu_ref)
                total_weights[model_ichans] +=  weights

            load_quiet = True
        for ipol in range(npol):
            aligned_port[ipol, np.where(total_weights > 0)[0]] /= \
                total_weights[np.where(total_weights > 0)[0]]
        model_port = aligned_port[0]
        niter -= 1
        count += 1
    if norm in ("mean", "max", "prof", "rms", "abs"):
        for ipol in range(npol):
            aligned_port[ipol] = normalize_portrait(aligned_port[ipol], norm,
                                                    weights=None)
    if rot_phase:
        aligned_port = rotate_data(aligned_port, rot_phase)
    if place is not None:
        prof = np.average(aligned_port[0], axis=0, weights=None)
        delta = prof.max() * gaussian_profile(len(prof), place, 0.0001)
        phase = fit_phase_shift(prof, delta, Ns=nbin).phase
        aligned_port = rotate_data(aligned_port, phase)
    arch = model_data.arch
    arch.tscrunch()
    if pscrunch:
        arch.pscrunch()
    else:
        arch.convert_state("Stokes")
    arch.set_dispersion_measure(0.0)
    for subint in arch:
        for ipol in range(arch.get_npol()):
            for ichan in range(arch.get_nchan()):
                prof = subint.get_Profile(ipol, ichan)
                prof.get_amps()[:] = aligned_port[ipol, ichan]
                if total_weights[ichan].sum() == 0.0:
                    subint.set_weight(ichan, 0.0)
                else:
                    # subint.set_weight(ichan, weight)
                    subint.set_weight(ichan, 1.0)
    arch.unload(outfile)
    if not quiet: print("\nUnloaded %s.\n" % outfile)


if __name__ == "__main__":

    from optparse import OptionParser

    usage = "Usage: %prog -M <metafile> [options]"
    parser = OptionParser(usage)
    # parser.add_option("-h", "--help",
    #                  action="store_true", dest="help", default=False,
    #                  help="Show this help message and exit.")
    parser.add_option("-M", "--metafile_band3",
                      default="band3.200.T.meta",
                      action="store", metavar="metafile_band3", dest="metafile_band3",
                      help="Metafile of archives to average together.")
    parser.add_option("-L", "--metafile_band5",
                      default="band5.200.T.meta",
                      action="store", metavar="metafile_band5", dest="metafile_band5",
                      help="Metafile of archives to average together.")
    parser.add_option("-I", "--init_band3",
                      default="J1939+2134_59540.597823_500.rfiClean.fits.pazi.rot",
                      action="store", metavar="initial_guess_band3",
                      dest="initial_guess_band3",
                      help="Archive containing initial alignment guess.  psradd is used if -I is not used.")
    parser.add_option("-J", "--init_band5",
                      default="J1939+2134_59540.597681_1460.rfiClean.fits.pazi.rot",
                      action="store", metavar="initial_guess_band5",
                      dest="initial_guess_band5",
                      help="Archive containing initial alignment guess.  psradd is used if -I is not used.")
    parser.add_option("-g", "--width",
                      default=None,
                      action="store", metavar="fwhm", dest="fwhm",
                      help="Use a single Gaussian component of given FWHM to align archives.  Overides -I.")
    parser.add_option("-D", "--no_DM",
                      default=True,
                      action="store_false", dest="fit_dm",
                      help="Align the subintegrations/archives with a fit for phase only.")
    parser.add_option("-T", "--tscr",
                      default=False,
                      action="store_true", dest="tscrunch",
                      help="Tscrunch archives for the iterations.  Recommended unless there is reason to keep subint resolution (may speed things up).")
    parser.add_option("-p", "--poln",
                      default=True,
                      action="store_false", dest="pscrunch",
                      help="Output average Stokes portraits, not just total intensity.  Archives are internally converted or skipped (if state == 'Intensity').")
    parser.add_option("-C", "--cutoff",
                      default=0.0,
                      action="store", metavar="SNR_cutoff", dest="SNR_cutoff",
                      help="S/N ratio cutoff to apply to input archives. [default=0.0]")
    parser.add_option("-o", "--outfile_band3",
                      default="J1939_band3.port",
                      action="store", metavar="outfile_band3", dest="outfile_band3",
                      help="Name of averaged output archive. [default=metafile.algnd.fits]")
    parser.add_option("-z", "--outfile_band5",
                      default="J1939_band5.port",
                      action="store", metavar="outfile_band5", dest="outfile_band5",
                      help="Name of averaged output archive. [default=metafile.algnd.fits]")
    parser.add_option("-P", "--palign",
                      default=False,
                      action="store_true", dest="palign",
                      help="Passes -P to psradd if -I is not used. [default=False]")
    parser.add_option("-N", "--norm",
                      action="store", metavar="normalization", dest="norm",
                      default=None,
                      help="Normalize the final averaged data by channel ('None' [default], 'mean', 'max' (not recommended), 'prof', 'rms', or 'abs').")
    parser.add_option("-s", "--smooth",
                      default=False,
                      action="store_true", dest="smooth",
                      help="Output a second averaged archive, smoothed with psrsmooth -W. [default=False]")
    parser.add_option("-r", "--rot",
                      default=0.0,
                      action="store", metavar="phase", dest="rot_phase",
                      help="Additional rotation to add to averaged archive. [default=0.0]")
    parser.add_option("--place",
                      default=None,
                      action="store", metavar="place", dest="place",
                      help="Roughly place pulse to be at the phase given.  Overrides --rot. [default=None]")
    parser.add_option("--niter",
                      action="store", metavar="int", dest="niter", default=1,
                      help="Number of iterations to complete. [default=1]")
    parser.add_option("--verbose",
                      action="store_false", dest="quiet", default=True,
                      help="More to stdout.")

    (options, args) = parser.parse_args()

    if options.metafile_band3 is None or options.metafile_band5 is None or not options.niter:
        print("\nppalign.py - Aligns and averages homogeneous archives by fitting DMs and phases\n")
        parser.print_help()
        print("")
        parser.exit()

    metafile_band3 = options.metafile_band3
    metafile_band5 = options.metafile_band5
    initial_guess_band3 = options.initial_guess_band3
    initial_guess_band5 = options.initial_guess_band5
    fwhm = options.fwhm
    fit_dm = options.fit_dm
    tscrunch = options.tscrunch
    pscrunch = options.pscrunch
    SNR_cutoff = float(options.SNR_cutoff)
    outfile_band3 = options.outfile_band3
    outfile_band5 = options.outfile_band5
    palign = options.palign
    norm = options.norm
    smooth = options.smooth
    rot_phase = np.float64(options.rot_phase)
    if options.place is not None:
        rot_phase = 0.0
        place = np.float64(options.place)
    else:
        place = None
    niter = int(options.niter)
    quiet = options.quiet

    rm = False
    if initial_guess_band3 is None and initial_guess_band5 is None and fwhm is None:  # use psradd :(
        tmp_file = "ppalign.%d.tmp.fits" % np.random.randint(100, 1000)
        psradd_archives(metafile_band3, outfile=tmp_file, palign=palign)
        psradd_archives(metafile_band5, outfile=tmp_file, palign=palign)
        initial_guess = tmp_file
        rm = True
    elif fwhm:  # use fixed Gaussian component
        tmp_file = "ppalign.%d.tmp.fits" % np.random.randint(100, 1000)
        archive_band3 = open(metafile_band3, 'r').readlines()[0][:-1]  # use first archive
        archive_band5 = open(metafile_band5, 'r').readlines()[0][:-1]  # use first archive
        vap_cmd_band3 = "vap -n -c nbin %s" % archive_band3
        vap_cmd_band5 = "vap -n -c nbin %s" % archive_band5
        nbin_band3 = int(sub.Popen(shlex.split(vap_cmd_band3), stdout=sub.PIPE
                             ).stdout.readlines()[0].split()[-1])
        nbin_band5 = int(sub.Popen(shlex.split(vap_cmd_band5), stdout=sub.PIPE
                             ).stdout.readlines()[0].split()[-1])
        profile_band3 = gaussian_profile(nbin_band3, 0.5, float(fwhm))
        profile_band5 = gaussian_profile(nbin_band5, 0.5, float(fwhm))
        make_constant_portrait(archive_band3, tmp_file, profile=profile_band3, DM=0.0,
                               dmc=False, weights=None, quiet=quiet)
        make_constant_portrait(archive_band5, tmp_file, profile=profile_band5, DM=0.0,
                               dmc=False, weights=None, quiet=quiet)
        initial_guess = tmp_file
        rm = True
    else:# use initial archive provided
        print("yesyse")        
        vap_cmd_band3 = "vap -n -c nchan %s" % initial_guess_band3
        vap_cmd_band5 = "vap -n -c nchan %s" % initial_guess_band5
        nchan_band3 = int(sub.Popen(shlex.split(vap_cmd_band3), stdout=sub.PIPE
                              ).stdout.readlines()[0].split()[-1])
        nchan_band5 = int(sub.Popen(shlex.split(vap_cmd_band5), stdout=sub.PIPE
                              ).stdout.readlines()[0].split()[-1])
        
        if nchan_band3 == 1:  # used constant average profile
            tmp_file = "ppalign.%d.tmp.fits" % np.random.randint(100, 1000)
            archive_band3 = open(metafile_band3, 'r').readlines()[0][:-1]  # use first arch
            make_constant_portrait(archive_band3, tmp_file, profile=None, DM=0.0,
                                   dmc=False, weights=None, quiet=quiet)
            initial_guess = tmp_file
            rm = True
        if nchan_band5 == 1:  # used constant average profile
            tmp_file = "ppalign.%d.tmp.fits" % np.random.randint(100, 1000)
            archive_band5 = open(metafile_band5, 'r').readlines()[0][:-1]  # use first arch
            make_constant_portrait(archive_band5, tmp_file, profile=None, DM=0.0,
                                   dmc=False, weights=None, quiet=quiet)
            initial_guess = tmp_file
            rm = True
    
    if pscrunch:
        state = 'Intensity'
        npol = 1
    else:
        state = 'Stokes'
        npol = 4
    try:
        model_data_band3 = load_data(initial_guess_band3, state=state, dedisperse=False,
                            dededisperse=False, tscrunch=True, pscrunch=pscrunch,
                            fscrunch=False, rm_baseline=True, flux_prof=False,
                            refresh_arch=True, return_arch=True, quiet=quiet)
        model_data_band5 = load_data(initial_guess_band5, state=state, dedisperse=False,
                            dededisperse=False, tscrunch=True, pscrunch=pscrunch,
                            fscrunch=False, rm_baseline=True, flux_prof=False,
                            refresh_arch=True, return_arch=True, quiet=quiet)
    except IndexError:
        print("%s: has npol = 1; need npol == 4." % initial_guess_band3)
        print("%s: has npol = 1; need npol == 4." % initial_guess_band5)
        sys.exit()
    
    # print(model_data_band3)
    P_band3 = model_data_band3.Ps[0]
    P_band5 = model_data_band5.Ps[0]
    freqsx_band3 = model_data_band3.freqs[0]
    freqsx_band5 = model_data_band5.freqs[0]
    # print(freqsx_band3)
    # freqsx_band5 = model_data_band5.freqs[0, model_data_band5.ok_ichans[0]]
    # print(freqsx_band5)
    bedelay_band3 = model_data_band3.backend_delay
    bedelay_band5 = model_data_band5.backend_delay
    MJDs_band3 = np.array([model_data_band3.epochs[isub].in_days() \
                      for isub in range(1)], dtype=np.double)
    MJDs_band5 = np.array([model_data_band5.epochs[isub].in_days() \
                      for isub in range(1)], dtype=np.double)
    MJD_band3 = MJDs_band3[0]
    MJD_band5 = MJDs_band5[0]
    nu_mean_band3 = freqsx_band3.mean()
    nu_mean_band5 = freqsx_band5.mean()
    phase_offset_band3 = 0.0
    phase_offset_band5_p3 = old_div(((MJD_band3 - MJD_band5)*24.*3600 + (bedelay_band3 - bedelay_band5)),P_band3)
    phase_offset_band5_p5 = old_div(((MJD_band3 - MJD_band5)*24.*3600 + (bedelay_band3 - bedelay_band5)),P_band5)
    phase_offset_band5 = old_div(((MJD_band3 - MJD_band5)*24.*3600 + (bedelay_band3 - bedelay_band5)),P_band5)
    #print("P3_PO",phase_offset_band5_p3)
    #print("P5_PO",phase_offset_band5_p5)
    #print(phase_offset_band5)
    #print(P_band3)
    #print(P_band5)
    
    #print('%.40f'% MJD_band3)
    #print('%.40f'% MJD_band5)
    #print('%.40f'% bedelay_band3)
    #print('%.40f'% bedelay_band5)
    model_port_temp_band3 = (model_data_band3.masks * model_data_band3.subints)[0, 0]
    model_port_temp_band5 = (model_data_band5.masks * model_data_band5.subints)[0, 0]
    # print("phase_offset",phase_offset)
    
    #print(nu_mean_band3)
    #print(np.shape(model_port_temp_band3))
    #print(np.shape(model_port_temp_band5))
    #print(np.shape(freqsx_band3))
    #print(np.shape(freqsx_band5))
    #print(P_band3)
    #print(P_band5)
    DM = model_data_band3.DM
    rot_port_band3 = rotate_data(model_port_temp_band3, 0.0, DM, P_band3, freqsx_band3,
                            nu_mean_band3)  # why not nu_fit?
    rot_port_band5 = rotate_data(model_port_temp_band5, phase_offset_band5, DM, P_band5, freqsx_band5,
                            nu_mean_band3)  # why not nu_fit?
    
    # plt.plot(figsize=(20,20))
    # plt.imshow(model_port_temp_band3,extent=(0,1,0,1), aspect='auto')
    # plt.title("portx_34")
    # plt.xlabel('nbin')
    # plt.ylabel('nchan')
    # plt.show()
    # plt.plot(figsize=(20,20))
    # plt.imshow(model_port_temp_band5,extent=(0,1,0,1), aspect='auto')
    # plt.title("portx_34")
    # plt.xlabel('nbin')
    # plt.ylabel('nchan')
    # plt.show()
    # plt.plot(figsize=(20,20))
    # plt.imshow(rot_port_band3,extent=(0,1,0,1), aspect='auto')
    # plt.title("portx_34")
    # plt.xlabel('nbin')
    # plt.ylabel('nchan')
    # plt.show()
    # plt.plot(figsize=(20,20))
    # plt.imshow(rot_port_band5,extent=(0,1,0,1), aspect='auto')
    # plt.title("portx_34")
    # plt.xlabel('nbin')
    # plt.ylabel('nchan')
    # plt.show()
    # print(ok)
    
    align_archives(metafile_band3, initial_guess=initial_guess_band3, phase_offset = phase_offset_band3,
            nu_mean_band = nu_mean_band3, fit_dm=fit_dm,
            tscrunch=tscrunch, pscrunch=pscrunch, SNR_cutoff=SNR_cutoff,
            outfile=outfile_band3, norm=norm, rot_phase=rot_phase, place=place,
            niter=niter, quiet=quiet)
    align_archives(metafile_band5, initial_guess=initial_guess_band5, phase_offset = phase_offset_band5,
            nu_mean_band = nu_mean_band3, fit_dm=fit_dm,
            tscrunch=tscrunch, pscrunch=pscrunch, SNR_cutoff=SNR_cutoff,
            outfile=outfile_band5, norm=norm, rot_phase=rot_phase, place=place,
            niter=niter, quiet=quiet)

    if smooth:
        if outfile is None:
            outfile = metafile + ".algnd.fits"
        psrsmooth_archive(outfile, options="-W")
    if rm:
        rm_cmd = "rm -f %s" % tmp_file
        rm_call = sub.Popen(shlex.split(rm_cmd))
        rm_call.wait()
