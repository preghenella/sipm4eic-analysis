#include "analysis_utils.h"
#include "fine_analysis.h"

const int       frame_size      = 1024; // [clock cycles]
const double    reference_delay = 10.;  // [ns]
namespace au = analysis_utils;

void
fine_analysis_test(std::string dirname = "/Users/nrubini/Analysis/ePIC/sipm4eic-analysis/Data/decoded/1/", std::vector<std::string> filenames = au::all_filenames)
{

  /** 
   ** POPULATE FRAMED DATA
   **/
  
  au::framed_data_t framed_data;
  au::populate_framed_data(framed_data, dirname, filenames, frame_size);  
  
  /** 
   ** POST-PROCESSING ANALYSIS 
   **/
      
  std::cout << " --- post processing " << std::endl;

  /** 
   ** HISTOGRAMS
   **/
      
  auto hN = new TH2F("hN", ";N_{TIME0} fired pixels;N_{TIME1} fired pixels;", 33, 0., 33., 33, 0., 33.);
  TH1 *hDelta[6] = {nullptr};
  TH2 *hMap[6] = {nullptr};
  for (int ichip = 0; ichip < 6; ++ichip) {
    hMap[ichip] = new TH2F(Form("hMap_%d", ichip), ";pixel row;pixel column", 8, 0., 8., 4, 0., 4.);
    hDelta[ichip] = new TH1F(Form("hDelta_%d", ichip), "hit - reference time (ns)", 2 * 60 * frame_size, -frame_size * au::coarse_to_ns, frame_size * au::coarse_to_ns);
  }
  
  /** 
   ** LOOP OVER DATA
   **/
      
  /** loop over spills **/
  for (auto &spill_data : framed_data) {
    auto spill = spill_data.first;
    auto &frames = spill_data.second;

    /** loop over frames **/
    for (auto &frame_data : frames) {
      auto frame = frame_data.first;
      auto &chips = frame_data.second;

      /** 
       ** COMPUTE REFERENCE TIME 
       **/
      
      double average_time[6] = {0.};

      /** loop over chips **/
      for (auto &chip_data : chips) {
        auto chip = chip_data.first;
        auto &channels = chip_data.second;
        
        /** loop over channels **/
        for (auto &channel_data : channels) {
          auto channel = channel_data.first;
          auto &hits = channel_data.second;

          /** compute time of first hit **/
          auto &hit = hits[0]; //( int fifo, int pixel, int column, int tdc, bool kUseFIFO );
            double time = hit.coarse * au::coarse_to_ns + hit.rollover * au::rollover_to_ns + calculate_calibrated_phase( hit.fine, "1", au::get_global_index( hit.fifo, hit.pixel, hit.column, hit.tdc ) ) * au::coarse_to_ns; // [ns]
          average_time[chip] += time;
          
        } /** end of loop over channels **/

        if (channels.size() > 0)
          average_time[chip] /= channels.size();
        
      } /** end of loop over chips **/
      
      /** request 32 fired pixels on both timing scintillators **/
      hN->Fill(chips[4].size(), chips[5].size());
      if (chips[4].size() != 32 || chips[5].size() != 32) continue;

      /** request compatible timing between timing scintillators **/      
      auto delta = average_time[4] - average_time[5];
      if (fabs(delta) > 5.) continue;

      /** compute reference time as average of scintillators 
          keeping in mind they are downstream and are delayed 
          with respect to cherenkov photons **/
      auto reference = 0.5 * (average_time[4] + average_time[5]) - reference_delay; // [ns]

      /** 
       ** CHERENKOV PHOTON TIMING 
       **/
      
      /** loop over chips **/
      for (auto &chip_data : chips) {
        auto chip = chip_data.first;
        auto &channels = chip_data.second;

        /** loop over channels **/
        for (auto &channel_data : channels) {
          auto channel = channel_data.first;
          auto &hits = channel_data.second;

          /** time of the first hit **/
          auto &hit = hits[0];
          auto index = au::get_index(hit.pixel, hit.column);
          double time = hit.coarse * au::coarse_to_ns + hit.rollover * au::rollover_to_ns; // [ns]
          double delta = time - reference;
          hDelta[chip]->Fill(delta);

          /** request hit time to be compatible with reference timing from scintillators **/
          if (std::fabs(delta) < 10.)
            hMap[chip]->Fill(index.first, index.second);
          
        } /** end of loop over channels **/

      } /** end of loop over chips **/
            
    } /** end of loop over frames **/
    
  } /** end of loop over spills **/

  /** 
   ** WRITE OUTPUT TO FILE
   **/
      
  auto fout = TFile::Open("analysis_example.root", "RECREATE");
  hN->Write();
  for (int ichip = 0; ichip < 6; ++ichip) {
    hMap[ichip]->Write();
    hDelta[ichip]->Write();
  }
  
  return;

}
