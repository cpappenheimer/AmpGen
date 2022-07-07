#include "AmpGen/Particle.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/QMI.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"



using namespace AmpGen;


real_t testSum(EventList_type& mcSig){
  real_t a =0;

  #pragma omp parallel for reduction(+:a)
  for(unsigned i=0;i<mcSig.size();++i) a+=1;
  return a; 
}





int main(int argc , char* argv[] ){
  OptionsParser::setArgs( argc, argv );
  MinuitParameterSet MPS;
  MPS.loadFromStream();

  TRandom3 rndm(0);

  EventType sigType(NamedParameter<std::string>("EventType", "", "Signal Type to generate"));
  size_t NInt(NamedParameter<size_t>("NInt", 1e7, "Number of events to calculate normalisation - should be large"));
  size_t seed(NamedParameter<size_t>("Seed", 0, "Random seed for generation"));
  size_t blockSize(NamedParameter<size_t>("blockSize", NInt, "Blocksize for generation"));
  size_t nEvents(NamedParameter<size_t>("nEvents", 1000, "number of events to generate, multiplies by the BR of each tag"));
  size_t nBins(NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms"));
  std::string outputFile(NamedParameter<std::string>("Output", "Generate_Output.root", "output file"));
  auto tags = NamedParameter<std::string>("TagTypes").getVector();
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  #ifdef _OPENMP
    omp_set_num_threads( nThreads );
    INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
    omp_set_dynamic( 0 );
  #endif




  auto  xy = QMI::xy(sigType, MPS); 
  EventList_type mcSig =  Generator<>(sigType, &rndm).generate(NInt);
  //real_t alpha(testSum(mcSig));
  //INFO("A = "<<alpha);
  //return 0;
  ProfileClock t_compilePoly, t_calcPoly, t_calcPolyManual;
  t_compilePoly.start();
  std::map<std::string, CompiledExpression<real_t(const real_t*, const real_t*)> > compiledPC(QMI::cPhaseCorrection(sigType, MPS));
  t_compilePoly.stop();
  INFO("Took "<<t_compilePoly.t_duration<<"ms to compile "<<NamedParameter<size_t>("PhaseCorrection::Order", 4)<<" polynomial");

  for (auto& p : compiledPC){
    INFO("Have "<<p.first<<"(0) = "<<p.second(mcSig[0].address())<<" x = "<<xy[0](mcSig[0].address())<<" y = "<<xy[1](mcSig[0].address()) ); 
  }
/*
  t_calcPoly.start();
//  real_t f0 = QMI::phaseCorrection(mcSig[0], compiledPC, MPS);
  real_t f0 =0;

  #if ENABLE_AVX  
        #pragma omp parallel for reduction(+:f0)
        for( unsigned i=0;i<mcSig.size();++i )//unsigned block = 0 ; block < m_events->nBlocks(); ++block )
        {
          //tmp[block] = m_events->weight(block) * AVX::log(this->operator()(m_events->block(block), block));  
          f0 += QMI::phaseCorrection(mcSig[i], compiledPC, MPS);
        }

  #endif
  


  t_calcPoly.stop();

  INFO("Took "<<t_calcPoly.t_duration<<"ms to calculate "<<f0);

  t_calcPolyManual.start();
  real_t f0_manual = 0;
//  real_t f0_manual = std::accumulate(
//    mcSig.begin(),
//    mcSig.end(),
//    0,
//    [&xy, &MPS](real_t rt, Event& evt){
//      return rt + QMI::phaseCorrection(evt, xy[0], xy[1], MPS);
//    }
//  );

  #if ENABLE_AVX 
        #pragma omp parallel for reduction(+:f0_manual)
        for( unsigned i=0;i<mcSig.size();++i)//unsigned block = 0 ; block < m_events->nBlocks(); ++block )
        {
          //tmp[block] = m_events->weight(block) * AVX::log(this->operator()(m_events->block(block), block));  
          f0_manual += QMI::phaseCorrection(mcSig[i], xy[0], xy[1], MPS);
        }

  #endif




  t_calcPolyManual.stop();

  INFO("Took "<<t_calcPolyManual.t_duration<<"ms to calculate "<<f0_manual);
*/
  

  MinuitParameterSet MPS_Kspipi;
  MPS_Kspipi.loadFromFile("Kspipi.opt");
  AddCPConjugate(MPS_Kspipi);
  CoherentSum A(sigType, MPS_Kspipi);

  CoherentSum Abar(sigType.conj(true), MPS_Kspipi);

  A.setEvents(mcSig);
  A.setMC(mcSig);
  A.prepare();
  Abar.setMC(mcSig);
  Abar.prepare();

  Event evt = mcSig[0];
  Event evtT = mcSig[0];
  evtT.swap(1,2);
  INFO("A(0) = "<<A.getValNoCache(evt));
  INFO("Abar(0) = "<<Abar.getValNoCache(evt));
  INFO("A(0T) = "<<A.getValNoCache(evtT));
  INFO("Abar(0T) = "<<Abar.getValNoCache(evtT));
  real_t dd = QMI::dd(evt, A, Abar);
  real_t ddT = QMI::dd(evtT, A, Abar);
  INFO("dd(0) = "<<dd);
  INFO("dd(0T) = "<<ddT);

  std::vector<Expression> Phi = QMI::dalitz(sigType);
  std::vector<CompiledExpression<real_t(const real_t*, const real_t*)> > cPhi;
  for (int i=0;i<Phi.size();i++){
      CompiledExpression<real_t(const real_t*, const real_t*)> cPhi_i(Phi[i], "Phi_" + std::to_string(i) , sigType.getEventFormat(), &MPS);
      cPhi_i.prepare(); cPhi_i.compile();
      cPhi.push_back(cPhi_i);
  }
  INFO("s23(0) = "<<cPhi[2](evt.address()));
  INFO("s23(0T) = "<<cPhi[2](evtT.address()));

  INFO("theta23(0) = "<<cPhi[3](evt.address()));
  INFO("theta23(0T) = "<<cPhi[3](evtT.address()));







  INFO("A(0) = "<<A.getVal(mcSig[0]));
  INFO("NA = "<<A.norm());
  std::map<std::string, std::vector<real_t> > amps(QMI::AmpArrays(mcSig, A, Abar));
/*
  Expression f = 1;
  Particle mother(sigType.decayDescriptor(), sigType.finalStates());
  Tensor p0(mother.P());
  Tensor p1(mother.daughter(0)->P());
  Tensor p2(mother.daughter(1)->P());
  Tensor p3(mother.daughter(2)->P());
  Expression s12 = 
    (p1[3] + p2[3]) * (p1[3] + p2[3]) - 
    (p1[0] + p2[0]) * (p1[0] + p2[0]) - 
    (p1[1] + p2[1]) * (p1[1] + p2[1]) - 
    (p1[2] + p2[2]) * (p1[2] + p2[2]);

 

    //using amp_type    = CompiledExpression<void(complex_v*, const size_t&, const real_t*, const float_v*)>;
    //CompiledExpression<complex_t(const real_t*, const real_t*)>(
      //    p.expression(), 
      //    p.decayDescriptor(),
      //    m_evtType.getEventFormat(), DebugSymbols(), disableBatch(), m_mps );

//  CompiledExpression<real_t(const real_t*, const float_v*)> comp;
  auto comp = CompiledExpression<real_t(const real_t*, const real_t*)> (s12, "s12", &MPS_Kspipi, sigType.getEventFormat());
  comp.prepare();
  comp.compile();
  real_t s = comp(mcSig[0].address()) ;
  INFO("s12(0) = "<<s);
*/

  TFile * tOutFile = TFile::Open(outputFile.c_str(), "RECREATE");
  tOutFile->cd();
  INFO("Calculate the cosine term for normalisation");
  real_t C(QMI::cosTerm(amps, MPS, compiledPC, mcSig));
  real_t S(QMI::sinTerm(amps, MPS, compiledPC, mcSig));
  auto my_A_real = [&A](Event& evt){ return std::real(A.getValNoCache(evt)); };
  auto my_A_imag = [&A](Event& evt){ return std::imag(A.getValNoCache(evt)); };
  auto my_Abar_real = [&Abar](Event& evt){ return std::real(Abar.getValNoCache(evt)); };
  auto my_Abar_imag = [&Abar](Event& evt){ return std::imag(Abar.getValNoCache(evt)); };


  for (auto tag : tags){
      INFO("tag = "<<tag);
      auto a = split(tag, ' ');



      EventType tagType(Particle(a[1], {}, false).eventType());
      EventList sigList(sigType);
      EventList tagList(tagType);
      EventList_type mcTag =  Generator<>(tagType, &rndm).generate(NInt);

      size_t nEvents_tag(std::stod(a[2]) * nEvents);

      bool sameType = sigType == tagType;
      if (!sameType){
        MinuitParameterSet MPS_Tag;
        INFO("Looking for "<<a[0] <<".opt");
        MPS_Tag.loadFromFile(a[0] + ".opt");
      
        CoherentSum B(tagType, MPS_Tag);

        CoherentSum Bbar(tagType.conj(true), MPS_Tag);

        B.setMC(mcTag);
        B.prepare();
        Bbar.setMC(mcTag);
        Bbar.prepare();

        std::map<std::string, std::vector<real_t> > BAmps (QMI::AmpArrays(mcTag, B, Bbar));
        real_t CB(QMI::cosTermNoCorrection(BAmps, mcTag));
        real_t SB(QMI::sinTermNoCorrection(BAmps, mcTag));

        
        //INFO("nB = "<<B.norm());
//        real_t norm(QMI::corr_norm(A, Abar, B, Bbar, MPS, mcSig, compiledPC, sameType));

        real_t norm = A.norm() * Bbar.norm() + Abar.norm() * B.norm() - 2 * C * CB - 2 * S * SB;
        //norm = norm * mcSig.size();

        INFO(a[0]<<" Norm term = "<<norm);
        //real_t probCorr_unnorm(Event& evt1, Event& evt2, CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, ce& x, ce& y, bool sameTag)
        auto psi_unnorm = [&A, &Abar, &B, &Bbar, &MPS, &compiledPC, &sameType, &norm](Event evt1, Event evt2){
          return QMI::probCorr_unnorm(evt1, evt2, A, Abar, B, Bbar, MPS, compiledPC, sameType)/norm;
        };
/*
        real_t norm_A_manual, norm_Abar_manual, norm_psi_manual;
        norm_A_manual =0;
        norm_Abar_manual =0;
        norm_psi_manual =0;
        #pragma omp parallel for reduction (+:norm_A_manual,norm_Abar_manual,norm_psi_manual)
        for (unsigned i=0;i<mcSig.size();++i) {
            norm_A_manual += std::norm(A.getValNoCache(mcSig[i]))/(real_t)mcSig.size();
            norm_Abar_manual += std::norm(Abar.getValNoCache(mcSig[i]))/(real_t)mcSig.size();
            norm_psi_manual += psi_unnorm(mcSig[i], mcTag[i]);
          
        }
        INFO("norm_A_manual = "<<norm_A_manual);
        INFO("A.norm() = "<<A.norm());
        INFO("norm_Abar_manual = "<<norm_Abar_manual);
        INFO("Abar.norm() = "<<Abar.norm());
        INFO("norm_psi_manual = "<<norm_psi_manual);
        INFO("norm(psi) = "<<norm);
*/
//        real_t max(getMax(mcSig))        

        QMI::generatePsi3770( psi_unnorm, sigList, tagList, nEvents_tag, seed, blockSize );
         



//        INFO("psi(0,0) = "<<psi_unnorm(mcSig[0], mcTag[0]));

      }
      else{
        //  auto psi_unnorm = [&A, &Abar, &A, &Abar, &MPS, &xy, &sameType](Event& evt1, Event& evt2)
        //    return QMI::probCorr_unnorm(evt1, evt2, A, Abar, A, Abar, MPS, xy[0], xy[1], sameType);
        //;


        //real_t corr_norm(CoherentSum& A, CoherentSum& Abar, CoherentSum& B, CoherentSum& Bbar, MinuitParameterSet& MPS, EventList& mcSig, std::map<std::string, ce>& cPij, bool sameTag)
        //real_t norm(QMI::corr_norm(A, Abar, A, Abar, MPS, mcSig, compiledPC, sameType));

    
        real_t norm =2 *  A.norm() * Abar.norm() - 2 * C * C - 2 * S * S;
        INFO(a[0]<<" Norm term = "<<norm);
        auto psi_unnorm = [&A, &Abar, &MPS, &compiledPC, &sameType, &norm](Event evt1, Event evt2){
          return QMI::probCorr_unnorm(evt1, evt2, A, Abar, A, Abar, MPS, compiledPC, sameType)/norm;
        };
        /*
        real_t norm_A_manual, norm_Abar_manual, norm_psi_manual;
        #pragma omp parallel for reduction (+:norm_A_manual,norm_Abar_manual,norm_psi_manual)
        for (unsigned i=0;i<mcSig.size();++i) {
            norm_A_manual += std::norm(A.getValNoCache(mcSig[i]))/(real_t)mcSig.size();
            norm_Abar_manual += std::norm(Abar.getValNoCache(mcSig[i]))/(real_t)mcSig.size();
            norm_psi_manual += psi_unnorm(mcSig[i], mcTag[i])/(real_t)mcSig.size();
            if (isnan(psi_unnorm(mcSig[i], mcTag[i]))) INFO("psi("<<i<<") is nan: ("<<mcSig[i].s(0, 1)<<", "<<mcSig[i].s(0, 2)<<"), ("<<mcTag[i].s(0, 1)<<", "<<mcTag[i].s(0, 2)<<")");
        }
        INFO("norm_A_manual = "<<norm_A_manual);
        INFO("A.norm() = "<<A.norm());
        INFO("norm_Abar_manual = "<<norm_Abar_manual);
        INFO("Abar.norm() = "<<Abar.norm());
        INFO("norm_psi_manual = "<<norm_psi_manual);
        INFO("norm(psi) = "<<norm);
        */

 
       QMI::generatePsi3770( psi_unnorm, sigList, tagList, nEvents_tag, seed, blockSize );
      }
//    for (auto& evt : sigList){
//      INFO("gauss("<<evt.s(0,1)<<", "<<evt.s(0, 2)<<" = "<<QMI::GaussDalitz(evt, MPS));
//      real_t A(MPS["PhaseCorrection::Gauss::A"]->mean());
//        real_t w_0(MPS["PhaseCorrection::Gauss::wminus_0"]->mean());
//        real_t muplus(MPS["PhaseCorrection::Gauss::mu_+"]->mean());
//        real_t sigmaplus(MPS["PhaseCorrection::Gauss::sigma_+"]->mean());
//        real_t muminus(MPS["PhaseCorrection::Gauss::mu_-"]->mean());
//        real_t sigmaminus(MPS["PhaseCorrection::Gauss::sigma_-"]->mean());
//        real_t s01(evt.s(0, 1));
//        real_t s02(evt.s(0, 2));
//        real_t c2 = -9.54231895051727e-05;
//        real_t m2 = 0.8051636393861085;
//        real_t z2 = (s01 - s02)/2;
//        real_t w2 = m2 * z2 + c2;
//
//
//        real_t gauss_pp = std::pow(s01 - muplus, 2);
//        real_t gauss_pm = std::pow(s01 - muminus, 2);
//        real_t gauss_mm = std::pow(s02 - muminus, 2);
//        real_t gauss_mp = std::pow(s02 - muplus, 2);
//        real_t exp_pp = gauss_pp/(2 * sigmaplus);
//        real_t exp_pm = gauss_pm/(2 * sigmaminus);
//        real_t exp_mp = gauss_mp/(2 * sigmaplus);
//        real_t exp_mm = gauss_mm/(2 * sigmaminus);
//        real_t g = 0;
//        if (s01>s02) g = A*std::erf(w2/w_0)*exp_pp * exp_mm;
//        if (s01<s02) g = A*std::erf(w2/w_0)*exp_pm * exp_mp;
//          INFO("A = "<<A<<" erf(w-/w0) = "<<std::erf(w2/w_0)
//          << " gausspp = "<<gauss_pp
//          << " exp_pp = "<<exp_pp
//          << " gausspm = "<<gauss_pm
//          << " exp_pm = "<<exp_pm
//          << " gaussmp = "<<gauss_mp
//          << " exp_mp = "<<exp_mp
//          << " gaussmm = "<<gauss_mm
//          << " exp_mm = "<<exp_mm
//          << " g = "<<g
//          );         
//      
//
//        }
//      

      
      
//      fillEventCustom(
//      sigList.tree(("Signal_" + a[0]).c_str())->Write();
//      tagList.tree(("Tag_" + a[0]).c_str())->Write();
//
//      auto sigPlots = sigList.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack));
//
//      for ( auto& plot : sigPlots ) plot->Write(("Signal_" + a[0] + "_" + plot->GetName()).c_str());
//      if( NamedParameter<bool>("plots_2d",true) == true ){
//        auto proj = sigList.eventType().defaultProjections(nBins);
//        for( size_t i = 0 ; i < proj.size(); ++i ){
//          for( size_t j = i+1 ; j < proj.size(); ++j ){ 
//            sigList.makeProjection( Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack), PlotOptions::Prefix(("Signal_" + a[0] + "_" + proj[i].name() + "_" + proj[j].name()).c_str())  )->Write( ); 
//          }
//        }
//      } 
//
//
//      auto tagPlots = tagList.makeDefaultProjections(PlotOptions::Bins(nBins), PlotOptions::LineColor(kBlack));
//      for ( auto& plot : tagPlots ) plot->Write(("Tag_" + a[0] + "_" + plot->GetName()).c_str());
//      if( NamedParameter<bool>("plots_2d",true) == true ){
//        auto proj = tagList.eventType().defaultProjections(nBins);
//        for( size_t i = 0 ; i < proj.size(); ++i ){
//          for( size_t j = i+1 ; j < proj.size(); ++j ){ 
//            tagList.makeProjection( Projection2D(proj[i], proj[j]), PlotOptions::LineColor(kBlack), PlotOptions::Prefix(("Tag_" + a[0] + "_" + proj[i].name() + "_" + proj[j].name()).c_str() ) )->Write( ); 
//          }
//        }
//      }
//      
//
//  
    if (NamedParameter<bool>("boostPsi3770", true)){
      EventList sigList2(sigList.eventType());
      EventList tagList2(tagList.eventType());
      QMI::boostPsi3770Events(sigList, tagList, sigList2, tagList2);
      QMI::writeEvents( sigList2, "Signal_" + a[0], nBins);
      QMI::writeEvents( tagList2, "Tag_" + a[0], nBins);
      QMI::writeEvents( sigList, "RF_Signal_" + a[0], nBins);
      QMI::writeEvents( tagList, "RF_Tag_" + a[0], nBins); 
    }
    else{

      QMI::writeEvents( sigList, "Signal_" + a[0], nBins);
      QMI::writeEvents( tagList, "Tag_" + a[0], nBins);
    }
    if (NamedParameter<bool>("writeAmpData", true)){
      QMI::writeValues(sigList, my_A_real, a[0] + "_A_real");
      QMI::writeValues(sigList, my_A_imag, a[0] + "_A_imag");
      QMI::writeValues(sigList, my_Abar_real, a[0] + "_Abar_real");
      QMI::writeValues(sigList, my_Abar_imag,a[0] + "_Abar_imag");
    }
      
  


  }
  auto my_dd = [&A, &Abar](Event& evt){
    return QMI::dd(evt, A, Abar);
  };
  auto my_corr = [&compiledPC, &MPS](Event& evt){
    return QMI::phaseCorrection(evt, compiledPC, MPS);
  };

  QMI::writeValues(mcSig, my_dd, "dd");
  QMI::writeValues(mcSig, my_corr, "corr");
  QMI::writeDalitz(mcSig);
  QMI::writeValues(mcSig, my_A_real, "A_real");
  QMI::writeValues(mcSig, my_A_imag, "A_imag");
  QMI::writeValues(mcSig, my_Abar_real, "Abar_real");
  QMI::writeValues(mcSig, my_Abar_imag, "Abar_imag");
  
  


  tOutFile->Close(); 



  return 0;
}
