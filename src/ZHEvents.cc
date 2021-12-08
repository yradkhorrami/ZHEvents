#include "ZHEvents.h"

using namespace lcio ;
using namespace marlin ;

ZHEvents aZHEvents;

ZHEvents::ZHEvents():
	Processor("ZHEvents"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	n_ZHDecays(0),
	n_ZDecays(0)
{
	_description = "ZHEvents checks the decay modes of Z/H bosons and accepts/rejects events wrt the decay mode and number of idolated leptons and jets";

	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"PfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollection",
					"Name of input jet collection",
					m_inputJetCollection,
					std::string("Durham_nJets")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"isoLepCollection",
					"Name of input Isolated Lepton collection",
					m_inputIsoLepCollection,
					std::string("ISOLeptons")
				);

	registerProcessorParameter(	"includ_bb",
					"Include Z/H->bb events",
					m_includ_bb,
					bool(true)
				);

	registerProcessorParameter(	"includ_cc",
					"Include Z/H->cc events",
					m_includ_cc,
					bool(true)
				);

	registerProcessorParameter(	"includ_ss",
					"Include Z/H->ss events",
					m_includ_ss,
					bool(true)
				);

	registerProcessorParameter(	"includ_uu",
					"Include Z/H->uu events",
					m_includ_uu,
					bool(true)
				);

	registerProcessorParameter(	"includ_dd",
					"Include Z/H->dd events",
					m_includ_dd,
					bool(true)
				);

	registerProcessorParameter(	"includ_gg",
					"Include Z/H->gg events",
					m_includ_gg,
					bool(true)
				);

	registerProcessorParameter(	"includ_ee",
					"Include Z/H->e+e- events",
					m_includ_ee,
					bool(true)
				);

	registerProcessorParameter(	"includ_mumu",
					"Include Z/H->mu+mu- events",
					m_includ_mumu,
					bool(true)
				);

	registerProcessorParameter(	"includ_tautau",
					"Include Z/H->tau+tau- events",
					m_includ_tautau,
					bool(true)
				);

	registerProcessorParameter(	"includ_nu1nu1",
					"Include Z/H->nu1nu1 events",
					m_includ_nu1nu1,
					bool(true)
				);

	registerProcessorParameter(	"includ_nu2nu2",
					"Include Z/H->nu2nu2 events",
					m_includ_nu2nu2,
					bool(true)
				);

	registerProcessorParameter(	"includ_nu3nu3",
					"Include Z/H->nu3nu3 events",
					m_includ_nu3nu3,
					bool(true)
				);

	registerProcessorParameter(	"includ_gammagamma",
					"Include Z/H->GammaGamma events",
					m_includ_gammagamma,
					bool(true)
				);

	registerProcessorParameter(	"includ_WW",
					"Include Z/H->W+W- events",
					m_includ_WW,
					bool(true)
				);

	registerProcessorParameter(	"includ_ZZ",
					"Include Z/H->ZZ events",
					m_includ_ZZ,
					bool(true)
				);

	registerProcessorParameter(	"includ_HH",
					"Include Z/H->HH events",
					m_includ_HH,
					bool(true)
				);

	registerProcessorParameter(	"includ_other",
					"Include Z/H->Other events",
					m_includ_other,
					bool(true)
				);

	registerProcessorParameter(	"includZee",
					"Include Z->e+e- events",
					m_includZe1e1,
					bool(true)
				);

	registerProcessorParameter(	"includZmumu",
					"Include Z->mu+mu- events",
					m_includZe2e2,
					bool(true)
				);

	registerProcessorParameter(	"includZtautau",
					"Include Z->tau+tau- events",
					m_includZe3e3,
					bool(true)
				);

	registerProcessorParameter(	"nJets",
					"Number of jet should be in the event",
					m_nJets,
					int(0)
				);

	registerProcessorParameter(	"nIsoLeps",
					"Number of Isolated Leptons should be in the event",
					m_nIsoLeps,
					int(0)
				);

	registerProcessorParameter(	"cheatDecayMode",
					"Cheat decay mode of boson fromMCTruth or use flavour-likeness / PFOType",
					m_cheatDecayMode,
					bool(true)
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("Output.root")
				);

}

void ZHEvents::init()
{

	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters();
	if ( m_fillRootTree )
	{
		m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
		m_pTTree = new TTree("SLDCorrection", "SLDCorrection");
		m_pTTree->SetDirectory(m_pTFile);
		h_ZHDecayMode = new TH1I( "ZHDecayMode" , "; Decay Mode" , 17 , -0.5 , 16.5 ); n_ZHDecays = 0;
		h_ZHDecayMode->GetXaxis()->SetBinLabel(1,"HH");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(2,"ZZ");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(3,"W^{+}W^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(4,"b#bar{b}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(5,"c#bar{c}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(6,"s#bar{s}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(7,"u#bar{u}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(8,"d#bar{d}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(9,"gg");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(10,"e^{+}e^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(11,"#mu^{+}#mu^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(12,"#tau^{+}#tau^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(13,"#nu_{e}#bar{#nu}_{e}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(14,"#nu_{#mu}#bar{#nu}_{#mu}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(15,"#nu_{#tau}#bar{#nu}_{#tau}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(16,"#gamma#gamma");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(17,"other");
		h_ZLeptonicDecayMode = new TH1I( "h_ZLeptonicDecayMode" , "; Decay Mode" , 3 , -0.5 , 2.5 ); n_ZDecays = 0;
		h_ZLeptonicDecayMode->GetXaxis()->SetBinLabel(1,"Z#rightarrow e^{+}e^{-}");
		h_ZLeptonicDecayMode->GetXaxis()->SetBinLabel(2,"Z#rightarrow #mu^{+}#mu^{-}");
		h_ZLeptonicDecayMode->GetXaxis()->SetBinLabel(3,"Z#rightarrow #tau^{+}#tau^{-}");
	}
	this->Clear();
}

void ZHEvents::Clear()
{

}

void ZHEvents::processRunHeader()
{
	streamlog_out(DEBUG) << "   processRunHeader called" << std::endl ;
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
	streamlog_out(DEBUG) << "   processRunHeader finished successfully" << std::endl ;
}

void ZHEvents::processEvent( EVENT::LCEvent *pLCEvent )
{
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	this->Clear();
	bool trueNJets = false;
	bool trueNIsoLeps = false;
	bool askedDecayMode = false;
	const EVENT::LCCollection *MCParticleCollection{};
	const EVENT::LCCollection *JetCollection{};
	const EVENT::LCCollection *IsoleptonCollection{};
	try
	{
		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		JetCollection = pLCEvent->getCollection( m_inputJetCollection );
		IsoleptonCollection = pLCEvent->getCollection( m_inputIsoLepCollection );

		int nJets = JetCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << nJets << " jets in event, looking for " << m_nJets << " jets" << std::endl;
		if ( nJets != m_nJets )
		{
			trueNJets = false;
			streamlog_out( DEBUG3 ) << "	Number of jets in the event mismatches the asked number of jets, --------EVENT REJECTED--------" << std::endl;
		}
		else
		{
			trueNJets = true;
			streamlog_out( DEBUG3 ) << "	Number of jets in the event matches the asked number of jets, --------EVENT ACCEPTED--------" << std::endl;
		}

		int nIsoLeps = IsoleptonCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << nIsoLeps << " issolated leptons in event, looking for " << m_nIsoLeps << " isolated leptons" << std::endl;
		if ( nIsoLeps != m_nIsoLeps )
		{
			trueNIsoLeps = false;
			streamlog_out( DEBUG3 ) << "	Number of isolated leptons in the event is less than the asked number of isolated leptons, --------EVENT REJECTED--------" << std::endl;
		}
		else
		{
			trueNIsoLeps = true;
			streamlog_out( DEBUG3 ) << "	Number of isolated leptons in the event matches the asked number of isolated leptons, --------EVENT ACCEPTED--------" << std::endl;
		}
		if ( m_cheatDecayMode )
		{
			int leptonicDecayMode = 0;
			int bosonDecayMode = 0;
			const EVENT::MCParticle *firstElectron = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( 4 ) );
			for ( int i_d = 0 ; i_d < firstElectron->getDaughters().size() ; ++i_d )
			{
				EVENT::MCParticle *daughter = firstElectron->getDaughters()[ i_d ];
				if ( abs( daughter->getPDG() ) == 23 || abs( daughter->getPDG() ) == 25 )
				{
					leptonicDecayMode = findDecayLeptonicMode( daughter );
				}
			}
		}
	}
	catch( DataNotAvailableException &e )
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

}

int ZHEvents::findDecayLeptonicMode( MCParticle *boson )
{
	int decayMode = -1;
	for ( int i_d1 = 0 ; i_d1 < boson->getDaughters().size() ; ++i_d1 )
	{
		EVENT::MCParticle *firstDaughter = boson->getDaughters()[ i_d1 ];
		for ( unsigned int i_type = 0 ; i_type < ( sizeof( leptonicPDGCodes ) / sizeof( *leptonicPDGCodes ) ) ; ++ i_type )
		{
			if ( firstDaughter->getPDG() == leptonicPDGCodes[ i_type ] )
			{
				for ( int i_d2 = 0 ; i_d2 < boson->getDaughters().size() ; ++i_d2 )
				{
					if ( i_d1 != i_d2 )
					{
						EVENT::MCParticle *secondDaughter = boson->getDaughters()[ i_d2 ];
						if ( firstDaughter->getPDG() == -1 * secondDaughter->getPDG() ) decayMode = i_type;
					}
				}
			}
		}
	}
	return decayMode;
}

int ZHEvents::findDecayMode( MCParticle *boson )
{
	for ( int i_d1 = 0 ; i_d1 < boson->getDaughters().size() ; ++i_d1 )
	{
		EVENT::MCParticle *firstDaughter = boson->getDaughters()[ i_d1 ];
		for ( unsigned int i_type = 0 ; i_type < ( sizeof( PDGCodes ) / sizeof( *PDGCodes ) ) ; ++ i_type )
		{
			if ( firstDaughter->getPDG() == PDGCodes[ i_type ] )
			{
				for ( int i_d2 = 0 ; i_d2 < boson->getDaughters().size() ; ++i_d2 )
				{
					if ( i_d1 != i_d2 )
					{
						EVENT::MCParticle *secondDaughter = boson->getDaughters()[ i_d2 ];
						if ( firstDaughter->getPDG() <= 16 && firstDaughter->getPDG() == -1 * secondDaughter->getPDG() )
						{
							return firstDaughter->getPDG();
						}
//						else if ( firstDaughter->getPDG() > 16 && firstDaughter->getPDG() == secondDaughter->getPDG() )
					}
				}
			}
		}
	}
}








int ZHEvents::isZHDecayedTo( const EVENT::LCCollection *MCParticleCollection , int parentPDG , int daughtersPDG , int &daughter1index , int &daughter2index )
{
	int elementFrom = 8;
	int elementTo = 20;
	int isDecaydToDaughter = 0;
	int d1PDG = abs( daughtersPDG );
	int d2PDG = ( abs( daughtersPDG ) < 20 ? -1 * d1PDG : d1PDG );
	for ( int i_d1 = elementFrom ; i_d1 < elementTo ; ++i_d1 )
	{
		const EVENT::MCParticle *daughter1 = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_d1 ) );
		if ( daughter1->getPDG() == d1PDG )
		{
			for ( unsigned int i_parent = 0 ; i_parent < daughter1->getParents().size() ; ++i_parent )
			{
				const EVENT::MCParticle *parent = daughter1->getParents()[ i_parent ];
				if ( parent->getPDG() == parentPDG )
				{
					for ( unsigned int i_d2 = 0 ; i_d2 < parent->getDaughters().size() ; ++i_d2 )
					{
						const EVENT::MCParticle *daughter2 = parent->getDaughters()[ i_d2 ];
						if ( daughter2 != daughter1 && daughter2->getPDG() == d2PDG )
						{
							isDecaydToDaughter = 1;
							daughter1index = i_d1;
							for ( int i_d = elementFrom ; i_d < elementTo ; ++i_d )
							{
								const EVENT::MCParticle *testMCP = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_d ) );
								if ( testMCP == daughter2 ) daughter2index = i_d;
							}
							return isDecaydToDaughter;
						}
					}
				}


				else if ( parent->getPDG() != 25 )
				{
					for ( unsigned int i_d2 = 0 ; i_d2 < parent->getDaughters().size() ; ++i_d2 )
					{
						const EVENT::MCParticle *daughter2 = parent->getDaughters()[ i_d2 ];
						if ( daughter2 != daughter1 && daughter2->getPDG() == d2PDG ) isDecaydToDaughter = 1;
					}
				}
			}
		}
	}
	return isDecaydToDaughter;
}

int ZHEvents::isZDecayedTo( const EVENT::LCCollection *MCParticleCollection , int parentPDG , int daughtersPDG , int daughter1index , int daughter2index )
{
	int elementFrom = 8;
	int elementTo = 20;
	int isDecaydToDaughter = 0;
	int d1PDG = abs( daughtersPDG );
	int d2PDG = ( abs( daughtersPDG ) < 20 ? -1 * d1PDG : d1PDG );
	for ( int i_d1 = elementFrom ; i_d1 < elementTo ; ++i_d1 )
	{
		const EVENT::MCParticle *daughter1 = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_d1 ) );
		if ( daughter1->getPDG() == d1PDG && i_d1 != daughter1index && i_d1 != daughter2index )
		{
			if ( daughter1->getParents().size() == 2 && abs( ( daughter1->getParents()[ 0 ] )->getPDG() ) == 11 && ( daughter1->getParents()[ 0 ] )->getPDG() == -1 * ( daughter1->getParents()[ 1 ] )->getPDG() )
			{
				const EVENT::MCParticle *parent = daughter1->getParents()[ 0 ];
				for ( unsigned int i_d2 = 0 ; i_d2 < parent->getDaughters().size() ; ++i_d2 )
				{
					const EVENT::MCParticle *daughter2 = parent->getDaughters()[ i_d2 ];
					if ( daughter2 != daughter1 && daughter2->getPDG() == d2PDG )
					{
						isDecaydToDaughter = 1;
						return isDecaydToDaughter;
					}
				}
			}
			else if ( daughter1->getParents().size() == 1 )
			{
				const EVENT::MCParticle *parent = daughter1->getParents()[ 0 ];
				if ( parent->getPDG() == parentPDG )
				{
					for ( unsigned int i_d2 = 0 ; i_d2 < parent->getDaughters().size() ; ++i_d2 )
					{
						const EVENT::MCParticle *daughter2 = parent->getDaughters()[ i_d2 ];
						if ( daughter2 != daughter1 && daughter2->getPDG() == d2PDG )
						{
							isDecaydToDaughter = 1;
							return isDecaydToDaughter;
						}
					}
				}
			}
		}
	}
	return isDecaydToDaughter;
}

void ZHEvents::check( EVENT::LCEvent *pLCEvent )
{
/*
	const EVENT::LCCollection *inJetCollection{};
	const EVENT::LCCollection *inIsoleptonCollection{};
	const EVENT::LCCollection *outPFOCollection{};
	const EVENT::LCCollection *outIsoleptonCollection{};
	try
	{
		inJetCollection = pLCEvent->getCollection( m_inputJetCollection );
		inIsoleptonCollection = pLCEvent->getCollection( m_inputIsoLepCollection );
		outPFOCollection = pLCEvent->getCollection( m_outputPfoCollection );
		outIsoleptonCollection = pLCEvent->getCollection( m_outputIsolepCollection );
		int nJets = inJetCollection->getNumberOfElements();
		int nInPFOs = 0;
		for ( int i_jet = 0 ; i_jet < nJets ; ++i_jet )
		{
			ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( inJetCollection->getElementAt( i_jet ) );
			nInPFOs += jet->getParticles().size();
		}
		int nInIsoLep = inIsoleptonCollection->getNumberOfElements();
		int nOutPFOs = outPFOCollection->getNumberOfElements();
		int nOutIsoLep = outIsoleptonCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << nJets << " jets with " << nInPFOs << " PFOs and " << nInIsoLep << " Isolated Leptons converted to " << nOutPFOs << " PFOs and " << nOutIsoLep << " Isolated Leptons" << std::endl;
	}
	catch( DataNotAvailableException &e )
        {
          streamlog_out( WARNING ) << "	Input/Output collections not found in event: " << m_nEvt << std::endl;
        }
*/
}

void ZHEvents::end()
{
	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		h_ZHDecayMode->Scale( 100.0 / n_ZHDecays );
		h_ZHDecayMode->Write();
		h_ZLeptonicDecayMode->Scale( 100.0 / n_ZDecays );
		h_ZLeptonicDecayMode->Write();
		m_pTFile->Close();
		delete m_pTFile;
	}
}
