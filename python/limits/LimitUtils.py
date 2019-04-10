#!/usr/bin/env python

import os, sys
import subprocess as sp
import glob
import argparse
import ROOT as rt

from SMSConfig import sms_models, VERSION

def do_command(cmd, no_exec=True):
    """
    Input: list of strings that should be joined to form the command.
    """
    print ' '.join(cmd)
    if not no_exec:
        sp.call(cmd)

def get_limit_dir(model, version=VERSION):
    return '/eos/cms/store/group/phys_susy/razor/Run2Analysis/Limits/RazorInclusive2016/{}/{}'.format(
            version, model)

def write_done_file(model, box, name='MADD'):
    """
    Creates text file containing filenames of completed jobs.
    Returns the name of the file.
    """
    done_jobs = glob.glob(get_limit_dir(model)
            +'/higgsCombine{}_*.root'.format(name))
    done_file_name = 'done_limits_{}_{}_{}.txt'.format(name, model, box)
    with open(done_file_name, 'w') as done_file:
        for f in done_jobs:
            done_file.write(os.path.basename(f)+'\n')
    return done_file_name

def submit_box_signif(model, boxes, 
        combined_with_boost=False, no_exec=True):
    script = 'python/ComputeSignificance.py'
    out_dir = get_limit_dir(model)
    queue = '1nh'
    combineName = 'MADDSignificance'
    if combined_with_boost:
        combineName = 'RazorInclusiveBoostSignficance'
    done_file_name = write_done_file(
            model, '_'.join(boxes), name=combineName)
    command = ['python', script, '--model', model,
            '--dir', out_dir, '--queue', queue,
            '--done-file', done_file_name, '--boxes']
    command += boxes
    if combined_with_boost:
        command.append('--combined-with-boost')
    if no_exec:
        command.append('--no-sub')
    do_command(command, False)

def submit_signif(model, tag, sms, combined=False, 
        combined_with_boost=False, no_exec=True):
    """
    Submit jobs to compute significance on existing cards
    """
    if combined or combined_with_boost:
        submit_box_signif(model, sms.boxes, combined_with_boost,
                no_exec)
    else:
        for box in sms.boxes:
            submit_box_signif(model, [box], False, no_exec)

def submit_box(model, box, bkg_dir=None, no_boost_cuts=False,
        fine_grained=True, no_sys=False, no_sub=True):
    """
    Submit jobs for a single box.
    Do not resubmit done jobs.
    """
    script = 'python/RunMADDLimitJobs.py'
    out_dir = VERSION
    queue = '8nh'
    done_file_name = write_done_file(model, box)
    command = ['python', script, '--box', box, '--model', model,
            '--dir', out_dir, '--queue', queue, 
            '--done-file', done_file_name]
    if bkg_dir is not None:
        command += ['--bkg-dir', bkg_dir]
    if no_sub:
        command.append('--no-sub')
    if no_boost_cuts:
        command.append('--no-boost-cuts')
    if fine_grained:
        command.append('--fine-grained')
    if no_sys:
        command.append('--no-sys')
    do_command(command, False)

def submit(model, tag, sms, bkg_dir=None, no_boost_cuts=False,
        fine_grained=False, no_sys=False, no_sub=True):
    """
    Submits the limit jobs for one SMS scan.
    model: string - name of model (without 'SMS-')
    tag: e.g. Razor2016_MoriondRereco
    sms: SMS object containing scan information
    no_boost_cuts: do not cut on number of t/w tags
    fine_grained: run on non-aggregated MC histograms
    no_sub: if True, do not submit jobs
    """
    for box in sms.boxes:
        print "Box {}".format(box)
        if sms.submodels is not None:
            for submodel in sms.submodels:
                print "Submitting jobs for {}".format(submodel)
                submit_box(submodel, box, bkg_dir, no_boost_cuts,
                        fine_grained, no_sys, no_sub)
        else:
            submit_box(model, box, bkg_dir, no_boost_cuts, 
                    fine_grained, no_sys, no_sub)

def submit_combine(model, tag, sms, no_sub=True):
    """
    Submits jobs to combine limits for one SMS scan.
    Combines all relevant boxes for the chosen model.
    """
    script = 'python/CombineMADDLimits.py'
    out_dir = get_limit_dir(model)
    queue = '1nd'
    done_file_name = write_done_file(model, '_'.join(sms.boxes))
    command = ['python', script, '--model', model,
            '--dir', out_dir, '--queue', queue,
            '--done-file', done_file_name, '--boxes']
    command += sms.boxes
    if no_sub:
        command.append('--no-sub')
    do_command(command, False)

def submit_boost_combine(model, tag, sms, no_sub=True):
    """
    Submits jobs to combine limits for one SMS scan.
    Combines razor inclusive limit with boost limit.
    """
    script = 'python/CombineInclusiveBoostLimits.py'
    out_dir = get_limit_dir(model)
    queue = '1nd'
    done_file_name = write_done_file(model, '_'.join(sms.boxes),
            name='RazorInclusiveBoost')
    command = ['python', script, '--model', model,
            '--dir', out_dir, '--queue', queue,
            '--done-file', done_file_name, '--boxes']
    command += sms.boxes
    if no_sub:
        command.append('--no-sub')
    do_command(command, False)

def aggregate(model, tag, sms, no_exec=True):
    """
    Makes a directory containing results from all submodels
    for a given SMS. Rename job files in the new directory
    to uniformize the output files.
    """
    if sms.submodels is None:
        print "No submodels to aggregate for {}".format(model)
        return
    out_dir = get_limit_dir(model)
    do_command(['mkdir', '-p', out_dir], no_exec)
    for submodel in sms.submodels:
        if submodel == model:
            print "Skipping submodel {}; name change unneeded".format(submodel)
            continue
        print "Dataset: {}".format(submodel)
        in_dir = get_limit_dir(submodel)
        in_files = glob.glob(in_dir
            +'/higgsCombineMADD_*.root')
        in_files += glob.glob(in_dir
            +'/RazorInclusiveMADD_*.*')
        for f in in_files:
            out_f = f.replace(submodel, model)
            if os.path.isfile(out_f):
                print "File {} already exists! Check for duplicate mass points."
            if f.endswith('.txt'):
                # change the name of the histogram file within the card
                sed_cmd = "sed 's/{}/{}/g' {} > {}".format(submodel, model,
                        f, out_f)
                print sed_cmd
                if not no_exec:
                    os.system(sed_cmd)
            else:
                do_command(['cp', f, out_f], no_exec)
        
def copy_lepton_limits_for_no_boost_cuts(model, no_exec=True):
    # This is supposed to work correctly whether or not
    # the current limit directory is the no-boost-cuts or boost-cuts one
    version = VERSION
    if version.endswith('_NoBoostCuts'):
        version = version.replace('_NoBoostCuts', '')
    boost_cuts_dir = get_limit_dir(model, version)
    no_boost_cuts_dir = get_limit_dir(model, version+'_NoBoostCuts')
    do_command(['mkdir', '-p', no_boost_cuts_dir], no_exec)
    lepton_boxes = ['LeptonMultiJet', 'LeptonSevenJet']
    for box in lepton_boxes:
        for f in glob.glob('{}/higgsCombineMADD_{}_SMS-{}_*.root'.format(
            boost_cuts_dir, box, model)):
            cmd = ['cp', f, no_boost_cuts_dir]
            do_command(cmd, no_exec)
        for f in glob.glob('{}/RazorInclusiveMADD_SMS-{}_*_*_{}.*'.format(
            boost_cuts_dir, model, box)):
            cmd = ['cp', f, no_boost_cuts_dir]
            do_command(cmd, no_exec)

def get_plot_dir():
    return 'PlotsSMS/plots/{}'.format(VERSION)

def get_config_dir():
    return 'PlotsSMS/config/{}'.format(VERSION)

def get(model, tag, sms, no_smooth=False, do_combined=False,
        do_boost_combined=False, no_exec=True, signif=False):
    """
    Retrieves limit results.  Args are the same as submit()
    except:
    no_smooth: do not fill gaps between points
    do_combined: retrieve combined limit 
    """
    get_script = 'python/GetCombineMADD.py'
    contour_script = 'python/Get2DContour.py'
    out_dir = get_limit_dir(model)
    boxes = sms.boxes
    if do_combined or do_boost_combined:
        boxes = ['_'.join(sms.boxes)]
    for box in boxes:
        print "Box {}".format(box)
        command = ['python', get_script, '--box', box, '--model', model]
        combine_name = None
        if do_boost_combined:
            combine_name = 'RazorInclusiveBoost'
            if signif:
                combine_name += 'Significance'
            if len(sms.boxes) == 1:
                combine_name += '_'+box
            comb_dir = out_dir+'/RazorInclusiveBoost'
            command += ['--in-dir', out_dir]
            command += ['--dir', comb_dir]
            do_command(['mkdir', '-p', comb_dir], no_exec)
        else:
            command += ['--dir', out_dir]
        if signif:
            command.append('--signif')
            if combine_name is None:
                combine_name = 'MADDSignificance'
        if combine_name is not None:
            command += ['--combine-name', combine_name]
        do_command(command, no_exec)
        if signif:
            return
        
        command = ['python', contour_script, '--box', box, '--model', model]
        if do_boost_combined:
            command += ['--dir', comb_dir]
        else:
            command += ['--dir', out_dir]
        if not sms.isGluino:
            command += ['--xsec-file', 'data/stop13TeV.txt']
        if no_smooth:
            command.append('--no-smooth')
        do_command(command, no_exec)
        if model == 'T2qq':
            # need extra 8x degenerate squark curve
            do_command(command+['--degen', '8'], no_exec)

def get_box_str(box):
    if 'Lepton' not in box:
        return 'razor_0L'
    for piece in box.split('_'):
        if piece == 'DiJet' or piece == 'MultiJet':
            return 'razor_0L+1L'
    return 'razor_1L'

def write_plot_config(model, box, blind=True, preliminary=True, 
        lumi=35900, no_exec=True, do_boost_combined=False):
    """
    Writes a config file for making SMS plot.
    Returns the name of the config.
    """
    in_dir = get_limit_dir(model)
    if do_boost_combined:
        in_dir += '/RazorInclusiveBoost'
    name = '{}_{}'.format(model, box)
    results_file = '{}/{}_results.root'.format(
            in_dir, name)
    if blind:
        hist_type = 'Exp'
    else:
        hist_type = 'Obs'
    hist_name = 'xsecUL_{}_{}'.format(hist_type, name)
    text = "HISTOGRAM {} {}\n".format(results_file, hist_name)
    text += "EXPECTED {} Exp_{} ExpPlus_{} ExpMinus_{} kRed kOrange\n".format(
            results_file, name, name, name)
    text += "EXPECTED2 {} Exp_{} ExpPlus2_{} ExpMinus2_{} kGray+1 kGray\n".format(
            results_file, name, name, name)
    if not blind:
        text += "OBSERVED {} Obs_{} ObsPlus_{} ObsMinus_{} kBlack kBlue-9\n".format(
                results_file, name, name, name)
        if model == 'T2qq':
            text += "OBSERVEDEXTRA {} Obs_{} ObsPlus_{} ObsMinus_{} kBlack kBlue-9\n".format(results_file.replace('.root', '_EXTRA.root'), name, name, name)
            text += "EXPECTEDEXTRA {} Exp_{} ExpPlus_{} ExpMinus_{} kRed kOrange\n".format(results_file.replace('.root', '_EXTRA.root'), name, name, name)
    text += "PRELIMINARY"
    if preliminary:
        text += " preliminary\n"
    else:
        text += " \n"
    text += "LUMI {}\n".format(lumi)
    text += "ENERGY 13\n"
    box_str = get_box_str(box)
    text += "BOXES {}\n".format(box_str)

    print "Plotting config:\n"
    print text
    out_dir = get_config_dir()
    out_f = '{}/{}.config'.format(out_dir, name)
    if blind:
        out_f = out_f.replace('.config', '_blinded.config')
    if preliminary:
        out_f = out_f.replace('.config', '_preliminary.config')
    if not no_exec:
        do_command(['mkdir', '-p', out_dir], no_exec)
        with open(out_f, 'w') as f:
            f.write(text)
        print "Plot config written to {}".format(out_f)
    return out_f

def plot(model, tag, sms, blind=True, preliminary=True, no_smooth=False,
        no_exec=True, do_combined=False, do_boost_combined=False):
    """
    Plots limit results.  Args are the same as submit() except for:
    blind: do not draw observed limit
    preliminary: write 'preliminary' on plots
    no_smooth: filename indicates that points were not smoothed.
    """
    plot_dir = get_plot_dir()
    do_command(['mkdir', '-p', plot_dir], no_exec)
    boxes = sms.boxes 
    if do_combined or do_boost_combined:
        boxes = ['_'.join(sms.boxes)]
    for box in boxes:
        print "Box: {}".format(box)
        config = write_plot_config(model, box, blind, preliminary, 
                no_exec=no_exec, do_boost_combined=do_boost_combined)
        name = '{}/{}_{}_'.format(plot_dir, model, box)
        if blind:
            name += 'Blinded_'
        if preliminary:
            name += 'Preliminary_'
        if no_smooth:
            name += 'NoSmooth_'
        if do_boost_combined:
            name += 'RazorBoostCombined_'
        script = 'PlotsSMS/python/makeSMSplots.py'
        command = ['python', script, config, name]
        do_command(command, no_exec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model')
    parser.add_argument('--tag', default='Razor2016_MoriondRereco')
    parser.add_argument('--no-exec', action='store_true')
    parser.add_argument('--no-sub', action='store_true',
            help='same as --no-exec')

    # actions
    parser.add_argument('--submit', action='store_true')
    parser.add_argument('--combined', action='store_true',
            help='do combination of boxes')
    parser.add_argument('--combined-with-boost', action='store_true',
            help='combine with razor boost analysis')
    parser.add_argument('--combined-hadronic-with-boost', action='store_true',
            help='combine hadronic box with boost analysis')
    parser.add_argument('--combined-hadronic', action='store_true',
            help='combine hadronic boxes only')
    parser.add_argument('--combined-leptonic', action='store_true',
            help='combine leptonic boxes only')
    parser.add_argument('--significance', action='store_true')
    parser.add_argument('--get', action='store_true')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--aggregate', action='store_true',
            help='combine datasets (for SMS split into multiple datasets')
    parser.add_argument('--finish', action='store_true',
            help='equivalent to --aggregate --get --plot')
    parser.add_argument('--copy-lepton-limits', action='store_true',
            help='copy lepton limit results to NoBoostCut folder')

    # configuration options
    parser.add_argument('--blind', action='store_true')
    parser.add_argument('--preliminary', action='store_true')
    parser.add_argument('--no-smooth', action='store_true',
            help='do not fill gaps between plotted points')
    parser.add_argument('--bkg-dir', help='Specify background histogram location')
    parser.add_argument('--no-boost-cuts', action='store_true')
    parser.add_argument('--no-fine-grained', action='store_true')
    parser.add_argument('--no-sys', action='store_true')

    args = parser.parse_args()
    fine_grained = not args.no_fine_grained
    no_exec = (args.no_exec or args.no_sub)

    args.preliminary = True

    print "Model: {}".format(args.model)
    print "Analysis tag: {}".format(args.tag)
    try:
        sms = sms_models[args.model]
    except KeyError:
        sys.exit("Model {} is not implemented!".format(args.model))

    if args.combined_hadronic:
        print "Combining hadronic boxes"
        sms.boxes = [b for b in sms.boxes 
                if b in ['DiJet', 'MultiJet', 'SevenJet']]
        args.combined = True
    elif args.combined_leptonic:
        print "Combining leptonic boxes"
        sms.boxes = [b for b in sms.boxes 
                if b in ['LeptonMultiJet', 'LeptonSevenJet']]
        args.combined = True
    elif args.combined_hadronic_with_boost:
        print "Combining boost analysis with hadronic box"
        sms.boxes = ['MultiJet', 'SevenJet']
        args.combined_with_boost = True

    if args.significance and args.submit:
        print "Compute significance for {}\n".format(args.model)
        submit_signif(args.model, args.tag, sms, 
                args.combined, args.combined_with_boost, 
                no_exec)
        sys.exit()
    if args.submit:
        print "Submit limit jobs for {}\n".format(args.model)
        if args.combined_with_boost:
            submit_boost_combine(args.model, args.tag, sms, no_exec)
        elif args.combined:
            submit_combine(args.model, args.tag, sms, no_exec)
        else:
            submit(args.model, args.tag, sms, args.bkg_dir, 
                    args.no_boost_cuts, fine_grained, args.no_sys, no_exec)
        sys.exit()

    if (args.aggregate or args.finish) and not (args.combined
            or args.combined_with_boost):
        if sms.submodels is not None:
            print "Combine limit jobs from all datasets for {}".format(
                    args.model)
            aggregate(args.model, args.tag, sms, no_exec)
        else:
            print "Ignoring --aggregate, no submodels to aggregate"

    if args.copy_lepton_limits:
        print "Copy lepton limit files to NoBoostCut folder"
        copy_lepton_limits_for_no_boost_cuts(args.model, no_exec)

    if args.get or args.finish:
        print "Get limit results for {}\n".format(args.model)
        get(args.model, args.tag, sms, no_smooth=args.no_smooth, 
                do_combined=args.combined, no_exec=no_exec,
                do_boost_combined=args.combined_with_boost,
                signif=args.significance)

    if args.plot or args.finish:
        print "Plot limit for {}\n".format(args.model)
        plot(args.model, args.tag, sms, blind=args.blind,
                preliminary=args.preliminary, no_smooth=args.no_smooth,
                do_combined=args.combined, no_exec=no_exec,
                do_boost_combined=args.combined_with_boost)
