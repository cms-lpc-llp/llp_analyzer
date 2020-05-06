#include "SimpleTable.h"
#include <Riostream.h>
#include <TSystem.h>
#include <TFormula.h>

//ClassImp(SimpleTable)
//ClassImp(SimpleTable::MyParameter)

//--------------------------------------------------------------------------------------------------
void SimpleTable::MyParameter::Print(Option_t */*option*/) const
{
  // Print this parameter.

  printf("%s -> %e\n", GetName(), GetVal());
}

//--------------------------------------------------------------------------------------------------
SimpleTable::SimpleTable(const char *input) 
  : fTable(TCollection::kInitHashTableCapacity, 1)
{
  // Constructor.

  fTable.SetOwner(kTRUE);

  TString ifile(gSystem->ExpandPathName(input));

  if (ifile.IsNull())
    return;

  std::ifstream in(ifile);
  if (!in.good()) 
    return;

  Char_t dummy[1024];
  TString name;
  TString value;
  while(!in.eof()) {
    in >> name;
    if ((name.IsNull()) || name.BeginsWith("#")) {
      in.getline(dummy,1024);
      continue;
    }
    in >> value;
    in.getline(dummy,1024);

    if (value.IsNull()) {
      Error("SimpleTable", "Value corresponding to name %s is null.", name.Data());
      continue;
    }
    TFormula fval("formula",value);
    MyParameter *par = new MyParameter(name,fval.Eval(0));
    fTable.Add(par);
  }
}

//--------------------------------------------------------------------------------------------------
Double_t SimpleTable::Get(const char *name) const
{
  // Get value corresponding to name.

  const MyParameter *p = dynamic_cast<const MyParameter*>(fTable.FindObject(name));
  
  if(!p)
    Fatal("Get", "Could not get value for given name %s", name);

  return p->GetVal();
}

//--------------------------------------------------------------------------------------------------
Double_t SimpleTable::Has(const char *name) const
{
  // Return true if table contains given name.

  TObject *o = fTable.FindObject(name);
  return (o!=0);
}

//--------------------------------------------------------------------------------------------------
void SimpleTable::Print(Option_t *opt) const
{
  // Print content of table.
    std::cout << "Printing table (option " << opt << ")" << std::endl;
  TIter iter(fTable.MakeIterator());
  const MyParameter *p = dynamic_cast<const MyParameter*>(iter.Next());
  while (p) {
    p->Print();
    p = dynamic_cast<const MyParameter*>(iter.Next());
  }
}

