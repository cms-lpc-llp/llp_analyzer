{
TString rfitpath("$ROOFITSYS/include");
TString path = gSystem->GetIncludePath();
path += "-I. -I$ROOTSYS/src -I";
path += rfitpath;
gSystem->SetIncludePath(path.Data());
gSystem->Load("../python/lib/libRazorRun2.so");
}
