from icecube import dataio, dataclasses, icetray, MuonGun, recclasses
from icecube import spline_reco
from icecube.hdfwriter import I3SimHDFWriter
from icecube.hdfwriter import I3HDFWriter
import numpy as np
import argparse, os, time, sys
from I3Tray import *
import glob


def only_pframes(frame):
    if reco:
        if frame.Has('LLHFit_cleaned_linefit_1_5'):
            return True
    else:
        if frame.Has('MuonEffectiveArea'):
            return True
        
def book_files(indir, infile, outfile, reco):
    input_files = glob.glob(indir + '*.i3.*')
    print('{} i3 files' .format(len(input_files)))

    tray = I3Tray()
    tray.Add('I3Reader', FilenameList=input_files)
    tray.Add(only_pframes)
    if reco:
        print('writing using HDFWriter')
        tray.AddSegment(I3HDFWriter,
                Output = outfile,
                Keys = ['MCMuon','I3EventHeader', 'MuonEffectiveArea', 'linefit', 
                    'LLHFit_cleaned_linefit_1_5', 'LLHFit_cleaned_linefit_1_5FitParams'],
                SubEventStreams = ['fullevent'])
    else:
        print('writing using SimHDFWriter')
        tray.AddSegment(I3SimHDFWriter,
                Output = outfile,
                Keys = ['MCMuon','I3EventHeader', 'MuonEffectiveArea'])
    tray.AddModule('TrashCan','can')
    tray.Execute()
    tray.Finish()
    return



if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--indir', type=str, help='Directory to Read all i3 files')
    p.add_argument('--infile', type=str, default=None, help='option to pass one file to book')
    p.add_argument('--outfile', type=str, help='output file name')
    p.add_argument('--reco', type=bool, default=False, help='output file name')
    args = p.parse_args()
   
    indir = args.indir
    infile = args.infile
    outfile = args.outfile
    reco = args.reco

    book_files(indir, infile, outfile, reco)
