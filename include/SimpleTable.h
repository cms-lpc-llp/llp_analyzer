//--------------------------------------------------------------------------------------------------
//
// SimpleTable
//
// Simple wrapper around a THashTable. It can be used to access constants via a given name.
// Entries must be read from a plain text file.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef CMSANA_UTILS_SIMPLETABLE_H
#define CMSANA_UTILS_SIMPLETABLE_H

#include <THashTable.h>
#include <TParameter.h>

class SimpleTable : public TObject
{
    public:
        class MyParameter : public TParameter<double>
    {
        public:
            MyParameter() {}
            MyParameter(const char *n, Double_t v) : 
                TParameter<double>(n,v), fHash(TString(n).Hash()) {}

            ULong_t Hash()                       const { return fHash; }
            void    Print(Option_t */*option*/="") const;

        protected:
            ULong_t fHash; //stored hash value

            //ClassDef(MyParameter,1) // Extension to TParameter<double>
    };

    public:
        SimpleTable(const char *input);
        ~SimpleTable() {}

        Double_t         Get(const char *name)   const;
        Double_t         Has(const char *name)   const;
        void             Print(Option_t *opt="") const;

    protected:
        THashTable      fTable;

        //ClassDef(SimpleTable, 1) // Simple table to keep numbers
};
#endif
