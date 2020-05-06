#include <fstream>
#include <sstream>
#include <iterator>
#include <assert.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeFormula.h"
#include "SimpleTable.h"
#include "TKey.h"
#include "TDirectoryFile.h"

using namespace std;


//get list of files to open, add normalization branch to the tree in each file
int main(int argc, char* argv[]) {

    //parse input list to get names of ROOT files
    if(argc < 5){
        cerr << "usage SkimNtuple inputList.txt <outputDirectory> <outputfileLabel> <Skim Cut String>" << endl;
        return -1;
    }
    string inputList(argv[1]);
    
    string outputDir = argv[2];
    string outputfileLabel = argv[3];
    string SkimCutString = argv[4];

    ifstream filein(inputList.c_str());
    string curFilename;
    vector<string> inputLines;
    while(getline(filein, curFilename)){
        if(curFilename.at(0) != '#') inputLines.push_back(curFilename); //'#' denotes a comment
        else cout << "(Skipping commented line in input)" << endl;
    }

    //open each ROOT file and add the normalization branch
    for(auto& line : inputLines){
        //parse input -- input lines should be in the form datasetName fileName
        istringstream buf(line);
        istream_iterator<std::string> beg(buf), end;
        vector<std::string> inputs(beg, end);
        
        string fileName = inputs[0];

        //create output file
	string outputfilename = Form("%s/%s_%s.root", outputDir.c_str(), 
				     fileName.substr(fileName.find_last_of("/")+1, fileName.find_last_of(".")-fileName.find_last_of("/")-1).c_str(),
				     outputfileLabel.c_str());
	cout << "Output file: " << outputfilename << "\n";
        TFile *outputFile = new TFile(outputfilename.c_str(), "RECREATE");
	
        //loop over all TTrees in the file and add the weight branch to each of them
        TFile *inputFile = new TFile(fileName.c_str(), "READ");
        assert(inputFile);
        inputFile->cd();
        TIter nextkey(inputFile->GetListOfKeys());
        TKey *key;
        TKey *previous = NULL;
        string dirName = "";

        //if the first key is a TDirectoryFile, go inside it and skim there (temporary hack for cloning a single directory)
        /*TKey *firstkey = (TKey*)nextkey();
        string className = firstkey->GetClassName();
        if(className.compare("TDirectoryFile") == 0){
            TDirectoryFile* dir = (TDirectoryFile*)firstkey->ReadObj();
            dirName = dir->GetName(); 
            outputFile->mkdir(dirName.c_str());
            cout << "Entering directory " << dirName << endl;
            nextkey = TIter(dir->GetListOfKeys());
        }
        else { //reset it
            nextkey.Reset();
        }*/
        //end temporary hack

        while((key = (TKey*)nextkey())){
            string className = key->GetClassName();
            cout << "Getting key from file.  Class type: " << className << endl;
            //I haven't found a solution to copy arbitrary objects into the new file.
            //For now we only care about histograms and TTrees, so we 
            //handle those as special cases.
            if(className.find("TH1") != string::npos) {
                outputFile->cd();
                TH1F *outHist = (TH1F*)key->ReadObj();
                cout << "Copying histogram " << outHist->GetName() << " into output file" << endl;
                outHist->Write(outHist->GetName());
                inputFile->cd();
                continue;
            }
            if(className.compare("TTree") != 0){
                cout << "Skipping key (not a TTree)" << endl;
                outputFile->cd();
                TObject *outObj = key->ReadObj();
		cout << "Name: " << outObj->GetName() << " " << outObj->GetTitle() << "\n";
                outObj->Write(outObj->GetTitle());
                inputFile->cd();
                continue;
            }

            //if this key has the same name as the previous one, it's an unwanted cycle and we skip it
            if(previous != NULL && strcmp(key->GetName(), previous->GetName()) == 0)
            {
                continue;
            }
            previous = key;

            TTree *inputTree = (TTree*)key->ReadObj();
            cout << "Processing tree " << inputTree->GetName() << endl;

            //create new normalized tree
            outputFile->cd(dirName.c_str());
            TTree *outputTree = inputTree->CloneTree(0);  
            cout << "Events in the ntuple: " << inputTree->GetEntries() << endl;

	    TTreeFormula *formula = new TTreeFormula("SkimCutString", SkimCutString.c_str(), inputTree);
	    int EventsPassed = 0;

            //store the weights
            for (int n=0;n<inputTree->GetEntries();n++) { 
	      if (n%1000000==0) cout << "Processed Event " << n << "\n";
                inputTree->GetEntry(n);

		bool passSkim = false;
		
		passSkim = formula->EvalInstance();

		if (passSkim) {
		  EventsPassed++;
		  outputTree->Fill(); 
		}
            }

	    delete formula;
	    cout << "Skim Efficiency : " << EventsPassed << " / " << inputTree->GetEntries() 
		 << " = " << float(EventsPassed ) / float(inputTree->GetEntries()) 
		 << " \n";

            //save
            outputTree->Write();
            inputFile->cd();
        }
        inputFile->Close();
        cout << "Closing output file." << endl;

        outputFile->Close();
        delete outputFile;
    }
}
