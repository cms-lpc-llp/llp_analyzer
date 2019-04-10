class Marker(object):
    pass

if __name__ == '__main__':



    import inspect, os, sys
    from optparse import OptionParser
    topDir = os.path.abspath(os.path.dirname(inspect.getsourcefile(Marker)))

    # seting up root in lxplus and checking if CMSSW environment is set (needed for some python stuff)
    ROOTscript = "/afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00/slc4_ia32_gcc34/root/bin/thisroot.sh"
    if not 'ROOTSYS' in os.environ and os.path.exists(ROOTscript):
        if not 'CMSSW_BASE' in os.environ:
            print "echo \"ERROR: set CMSSW environment first\";\n"
            sys.exit();
        else:
            print "source " + ROOTscript + ";\n"
    elif not 'ROOTSYS' in os.environ:
        print "echo \"WARNING: could not setup root 5.28\";\n"

    parser = OptionParser()
    parser.add_option('-c',"--csh",dest="csh",action="store_true",default=False,
                  help="Use the csh rather than bash")
    options,args = parser.parse_args()

    def export_(var, string, csh):
        if csh:
            return 'setenv %s %s;' % (var,string)
        else:
            return 'export %s=%s;' % (var,string)
    def prepend_(var, string, csh):
        if var in os.environ:
            before = os.environ[var]
            if before.endswith(':'):
                before = before[0:-1]
            if csh:
                return 'setenv %s %s:%s;' % (var,string,before)
            else:
                return 'export %s=%s:${%s};' % (var,string,var)
        else:
            return export_(var,string,csh)
    print export_('RAZORFIT_BASE',topDir,options.csh),'\n'
    print prepend_('PYTHONPATH',os.path.join(topDir,'python'),options.csh),'\n'
    
    LIBDIR = os.path.join(topDir,'lib')
    #annoying difference between macos and linux...        
    if sys.platform.lower() in ['darwin']:
        print prepend_('DYLD_LIBRARY_PATH',LIBDIR,options.csh),'\n'
    else:
        print prepend_('LD_LIBRARY_PATH',LIBDIR,options.csh),'\n'
