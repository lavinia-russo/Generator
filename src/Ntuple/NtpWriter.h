//____________________________________________________________________________
/*!

\class   genie::NtpWriter

\brief   A simple class to facilitate creating the GENIE MC Ntuple from the
         output GENIE STDHEP event records.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 1, 2004

*/
//____________________________________________________________________________

#ifndef _NTP_WRITER_H_
#define _NTP_WRITER_H_

class TFile;
class TTree;

namespace genie {

class EventRecord;
class NtpMCEvent;

class NtpWriter {

public :

  NtpWriter();
  ~NtpWriter();
  
  void InitTree        (const char * filename);
  void AddStdhepRecord (int ievent, const EventRecord * ev_rec);
  void SaveTree        (void);

private:

  void Init(void);

  TFile *      fOutFile;
  TTree *      fOutTree;
  NtpMCEvent * fNtpMCEvent;
};

}      // genie namespace

#endif // _NTP_WRITER_H_
