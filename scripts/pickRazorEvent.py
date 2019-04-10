### Utility script to find out which big razor ntuple file contains a given event.
### Requires that you have a grid proxy set up

import sys
import argparse
import subprocess as sp
import glob
import ROOT as rt

def matches(name, dataset, version):
    """Returns true if we should use this dataset"""
    required = ['MINIAOD', dataset, version]
    for r in required:
        if r not in name:
            return False
    return True

def query_dataset_name(args):
    """
    Performs a DAS query to find the dataset 
    containing the target event
    """
    query = "'dataset run=%s'"%(args.run)
    cmd = 'das_client --query='+query+' --limit=0'
    print cmd
    try:
        output = sp.check_output(cmd, shell=True)
    except sp.CalledProcessError:
        print "Query didn't work.  Make sure you have a grid proxy set up"
        return
    fn = lambda name: matches(name, args.dataset, args.version)
    filtered = filter(fn, output.split('\n'))
    if len(filtered) == 0:
        print "No matches found!"
        return
    dataset = filtered[0]
    return dataset

def format_dataset_name(name):
    split = name.split('/')
    dataset = split[1]
    info = split[2].split('-')
    return '{}_{}_{}'.format(dataset, info[0].replace('Run',''), info[1])

def get_list_files(list_dir, dataset_name):
    return glob.glob('{}/{}*.cern.txt'.format(list_dir, dataset_name))

def found_in_file(fname, run, lumi, event):
    tfile = rt.TFile.Open(fname.replace('root://eoscms/', '').strip())
    tree = tfile.Get('ntuples/RazorEvents')
    tree.SetBranchStatus('*', 0)
    tree.SetBranchStatus('runNum', 1)
    tree.SetBranchStatus('lumiNum', 1)
    tree.SetBranchStatus('eventNum', 1)
    draw_str = 'runNum == {} && eventNum == {} && lumiNum == {}'.format(
            run, event, lumi)
    match = tree.Draw('>>elist1', draw_str)
    if match:
        pick_event(tree, run, lumi, event)
        return True
    return False

def pick_event(tree, run, lumi, event):
    elist = rt.gDirectory.Get("elist1")
    tree.SetBranchStatus("*", 1)
    tree.GetEntry(elist.GetEntry(0))
    out_f = rt.TFile.Open('razorevent_{}_{}_{}.root'.format(run, lumi, event),
            'RECREATE')
    out_tree = tree.GetTree().CloneTree(0)
    out_tree.Fill()
    out_tree.Write()
    out_f.Close()

def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('run')
    parser.add_argument('lumi')
    parser.add_argument('event')
    parser.add_argument('--dataset', default='HTMHT')
    parser.add_argument('--version', default='03Feb2017')
    return parser

if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()
    list_dir = 'lists/Run2/razorNtuplerV3p15/data/'
    
    dataset_name = query_dataset_name(args)
    if dataset_name:
        print "Event found in dataset {}".format(dataset_name)
    else:
        sys.exit("Event not found in DAS")
    dataset_list_name = format_dataset_name(dataset_name)
    for list_file in get_list_files(list_dir, dataset_list_name):
        print "Checking in list file {}".format(list_file)
        with open(list_file) as f:
            for line in f:
                print "Opening file {}".format(line)
                if found_in_file(line, args.run, args.lumi, args.event):
                    print "Event found in file: {}".format(line)
                    sys.exit()
    print "Sorry, event not found."
