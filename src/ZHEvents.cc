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
	m_nRecoIsoLeptons(0),
	m_nRecoJets(0),
	m_nTrueIsoLeptons(0),
	m_nTrueQuarks(0),
	m_leptonicDecayMode(0),
	m_bosonDecayMode(0),
	m_diLeptonTrueInvMass(0.0),
	m_diLeptonRecoInvMass(0.0),
	m_diQuarkInvMass(0.0),
	m_isoLepInvMassCutMin(0.0),
	m_isoLepInvMassCutMax(0.0),
	n_ZHDecays(0.0),
	n_ZDecays(0.0)
{
	_description = "ZHEvents checks the decay modes of Z/H bosons and accepts/rejects events wrt the decay mode and number of idolated leptons and jets";

	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
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

	registerInputCollection(	LCIO::MCPARTICLE,
					"trueQuarkCollection",
					"Name of input true quark collection",
					m_inputTrueQuarkCollection,
					std::string("trueQuarks")
				);

	registerInputCollection(	LCIO::MCPARTICLE,
					"trueIsoLepCollection",
					"Name of input true Isolated Lepton collection",
					m_inputTrueIsoLepCollection,
					std::string("TrueIsoLeps")
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

	registerProcessorParameter(	"isoLepInvMassCutMin",
					"Minimum Invariant mass of Isolated Leptons to accept events with diLeptonInvMass >= diLeptonInvMassCutMin",
					m_isoLepInvMassCutMin,
					float(0.0)
				);

	registerProcessorParameter(	"isoLepInvMassCutMax",
					"Maximum Invariant mass of Isolated Leptons to accept events with diLeptonInvMass <= diLeptonInvMassCutMax",
					m_isoLepInvMassCutMax,
					float(0.0)
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
		m_pTTree->Branch("nTrueIsoLeptons", &m_nTrueIsoLeptons , "nTrueIsoLeptons/I" );
		m_pTTree->Branch("nTrueQuarks", &m_nTrueQuarks , "nTrueQuarks/I" );
		m_pTTree->Branch("nRecoIsoLeptons", &m_nRecoIsoLeptons , "nRecoIsoLeptons/I" );
		m_pTTree->Branch("nRecoJets", &m_nRecoJets , "nRecoJets/I" );
		m_pTTree->Branch("leptonicDecayMode", &m_leptonicDecayMode , "leptonicDecayMode/I" );
		m_pTTree->Branch("bosonDecayMode", &m_bosonDecayMode , "bosonDecayMode/I" );
		m_pTTree->Branch("diLeptonTrueInvMass", &m_diLeptonTrueInvMass , "diLeptonTrueInvMass/F" );
		m_pTTree->Branch("diLeptonRecoInvMass", &m_diLeptonRecoInvMass , "diLeptonRecoInvMass/F" );
		m_pTTree->Branch("diQuarkInvMass", &m_diQuarkInvMass , "diQuarkInvMass/F" );
		h_ZHDecayMode = new TH1F( "ZHDecayMode" , "; Decay Mode" , 15 , -0.5 , 14.5 ); n_ZHDecays = 0.0;
		h_ZHDecayMode->GetXaxis()->SetBinLabel(1,"other");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(2,"W^{+}W^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(3,"b#bar{b}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(4,"c#bar{c}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(5,"s#bar{s}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(6,"u#bar{u}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(7,"d#bar{d}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(8,"gg");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(9,"e^{+}e^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(10,"#mu^{+}#mu^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(11,"#tau^{+}#tau^{-}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(12,"#nu_{e}#bar{#nu}_{e}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(13,"#nu_{#mu}#bar{#nu}_{#mu}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(14,"#nu_{#tau}#bar{#nu}_{#tau}");
		h_ZHDecayMode->GetXaxis()->SetBinLabel(15,"#gamma#gamma");
		h_ZLeptonicDecayMode = new TH1F( "ZLeptonicDecayMode" , "; Decay Mode" , 3 , 0.5 , 3.5 ); n_ZDecays = 0.0;
		h_ZLeptonicDecayMode->GetXaxis()->SetBinLabel(1,"Z#rightarrow e^{+}e^{-}");
		h_ZLeptonicDecayMode->GetXaxis()->SetBinLabel(2,"Z#rightarrow #mu^{+}#mu^{-}");
		h_ZLeptonicDecayMode->GetXaxis()->SetBinLabel(3,"Z#rightarrow #tau^{+}#tau^{-}");
	}
	this->Clear();
}

void ZHEvents::Clear()
{
	m_nTrueIsoLeptons = 0;
	m_nTrueQuarks = 0;
	m_nRecoIsoLeptons = 0;
	m_nRecoJets = 0;
	m_leptonicDecayMode = 0;
	m_bosonDecayMode = 0;
	m_diLeptonTrueInvMass = 0.0;
	m_diLeptonRecoInvMass = 0.0;
	m_diQuarkInvMass = 0.0;
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
	bool typeIsoLeps = false;
	bool trueIsoLepInvMass = false;
	bool recoIsoLepInvMass = false;
	bool askedDecayMode = false;
	bool askedLeptonicDecayMode = false;
	const EVENT::LCCollection *MCParticleCollection{};
	const EVENT::LCCollection *JetCollection{};
	const EVENT::LCCollection *IsoleptonCollection{};
	const EVENT::LCCollection *TrueQuarkCollection{};
	const EVENT::LCCollection *TrueIsoleptonCollection{};
	try
	{
		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		streamlog_out( DEBUG4 ) << "	MCParticle Collection: " << MCParticleCollection << std::endl;
		TrueQuarkCollection = pLCEvent->getCollection( m_inputTrueQuarkCollection );
		streamlog_out( DEBUG4 ) << "	True Quark Collection: " << TrueQuarkCollection << std::endl;
		TrueIsoleptonCollection = pLCEvent->getCollection( m_inputTrueIsoLepCollection );
		streamlog_out( DEBUG4 ) << "	True Isolated Lepton Collection: " << TrueIsoleptonCollection << std::endl;
		JetCollection = pLCEvent->getCollection( m_inputJetCollection );
		streamlog_out( DEBUG4 ) << "	Jet Collection: " << JetCollection << std::endl;
		IsoleptonCollection = pLCEvent->getCollection( m_inputIsoLepCollection );
		streamlog_out( DEBUG4 ) << "	Isolated Lepton Collection: " << IsoleptonCollection << std::endl;

		m_nRecoJets = JetCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << m_nRecoJets << " jets in event, looking for " << m_nJets << " jets" << std::endl;
		if ( m_nRecoJets != m_nJets )
		{
			trueNJets = false;
			streamlog_out( DEBUG3 ) << "	Number of jets in the event mismatches the asked number of jets, --------EVENT REJECTED--------" << std::endl;
		}
		else
		{
			trueNJets = true;
			streamlog_out( DEBUG3 ) << "	Number of jets in the event matches the asked number of jets, --------EVENT ACCEPTED--------" << std::endl;
		}
		setReturnValue( "trueNJets" , trueNJets );

		m_nRecoIsoLeptons = IsoleptonCollection->getNumberOfElements();
		streamlog_out( DEBUG4 ) << "	" << m_nRecoIsoLeptons << " isolated leptons in event, looking for " << m_nIsoLeps << " isolated leptons" << std::endl;
		if ( m_nRecoIsoLeptons != m_nIsoLeps )
		{
			trueNIsoLeps = false;
			streamlog_out( DEBUG3 ) << "	Number of isolated leptons in the event is less than the asked number of isolated leptons, --------EVENT REJECTED--------" << std::endl;
		}
		else
		{
			trueNIsoLeps = true;
			streamlog_out( DEBUG3 ) << "	Number of isolated leptons in the event matches the asked number of isolated leptons, --------EVENT ACCEPTED--------" << std::endl;
		}
		int isoLepType = 0;
		TLorentzVector recoDiLepton4Momentum( 0.0 , 0.0 , 0.0 , 0.0 );
		for ( int i_lep = 0 ; i_lep < m_nRecoIsoLeptons ; ++i_lep )
		{
			ReconstructedParticle *lepton = dynamic_cast<ReconstructedParticle*>( IsoleptonCollection->getElementAt( i_lep ) );
			recoDiLepton4Momentum += TLorentzVector( lepton->getMomentum()[ 0 ] , lepton->getMomentum()[ 1 ] , lepton->getMomentum()[ 2 ] , lepton->getEnergy() );
			isoLepType += lepton->getType();
		}
		m_diLeptonRecoInvMass = recoDiLepton4Momentum.M();
		recoIsoLepInvMass = ( recoDiLepton4Momentum.M() >= m_isoLepInvMassCutMin && recoDiLepton4Momentum.M() <= m_isoLepInvMassCutMax ? true : false );
		typeIsoLeps = ( isoLepType == 0 ? true : false );
		trueNIsoLeps = trueNIsoLeps && typeIsoLeps;
		setReturnValue( "trueNIsoLeps" , trueNIsoLeps );
		setReturnValue( "recoIsoLepInvMass" , recoIsoLepInvMass );

		m_nTrueIsoLeptons = TrueIsoleptonCollection->getNumberOfElements();
		TLorentzVector trueDiLepton4Momentum( 0.0 , 0.0 , 0.0 , 0.0 );
		for ( int i_lep = 0 ; i_lep < m_nTrueIsoLeptons ; ++i_lep )
		{
			MCParticle* lepton = dynamic_cast<MCParticle*>( TrueIsoleptonCollection->getElementAt( i_lep ) );
			trueDiLepton4Momentum += TLorentzVector( lepton->getMomentum()[ 0 ] , lepton->getMomentum()[ 1 ] , lepton->getMomentum()[ 2 ] , lepton->getEnergy() );
		}
		if ( m_nTrueIsoLeptons == m_nIsoLeps ) m_diLeptonTrueInvMass = trueDiLepton4Momentum.M();
		trueIsoLepInvMass = ( trueDiLepton4Momentum.M() >= m_isoLepInvMassCutMin && trueDiLepton4Momentum.M() <= m_isoLepInvMassCutMax ? true : false );
		setReturnValue( "trueIsoLepInvMass" , trueIsoLepInvMass );

		m_nTrueQuarks = TrueQuarkCollection->getNumberOfElements();
		TLorentzVector diQuark4Momentum( 0.0 , 0.0 , 0.0 , 0.0 );
		for ( int i_quark = 0 ; i_quark < m_nTrueQuarks ; ++i_quark )
		{
			MCParticle* quark = dynamic_cast<MCParticle*>( TrueQuarkCollection->getElementAt( i_quark ) );
			diQuark4Momentum += TLorentzVector( quark->getMomentum()[ 0 ] , quark->getMomentum()[ 1 ] , quark->getMomentum()[ 2 ] , quark->getEnergy() );
		}
		if ( m_nTrueQuarks == m_nJets ) m_diQuarkInvMass = diQuark4Momentum.M();

		int leptonicDecayMode = 0;
		int bosonDecayMode = 0;
		if ( m_cheatDecayMode )
		{
			const EVENT::MCParticle *leptonicZ = NULL;
			mcpVector leptonicDecayProducts{};
			mcpVector otherDecayProducts{};
			const EVENT::MCParticle *firstElectron = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( 4 ) );
			leptonicDecayMode = findDecayLeptonicMode( firstElectron , leptonicDecayProducts );
			if ( leptonicDecayMode <= 0 )
			{
				for ( unsigned int i_d = 0 ; i_d < firstElectron->getDaughters().size() ; ++i_d )
				{
					const EVENT::MCParticle *daughter = firstElectron->getDaughters()[ i_d ];
					if ( abs( daughter->getPDG() ) == 23 )
					{
						leptonicDecayMode = findDecayLeptonicMode( daughter , leptonicDecayProducts );
						if ( leptonicDecayMode > 0 ) leptonicZ = daughter;
					}
				}
				leptonicDecayProducts.clear();
			}
			if ( leptonicZ != NULL ) leptonicDecayMode = findDecayLeptonicMode( leptonicZ , leptonicDecayProducts );
			if ( m_nTrueIsoLeptons == 2 )
			{
				MCParticle* lepton = dynamic_cast<MCParticle*>( TrueIsoleptonCollection->getElementAt( 0 ) );
				int leptonPDG = abs( lepton->getPDG() );
				if ( leptonPDG == 11 ) leptonicDecayMode = 1;
				if ( leptonPDG == 13 ) leptonicDecayMode = 2;
				if ( leptonPDG == 15 ) leptonicDecayMode = 3;
			}
			if ( m_includZe1e1 && leptonicDecayMode == 1 ) askedLeptonicDecayMode = true;
			else if ( m_includZe2e2 && leptonicDecayMode == 2 ) askedLeptonicDecayMode = true;
			else if ( m_includZe3e3 && leptonicDecayMode == 3 ) askedLeptonicDecayMode = true;
			else  askedLeptonicDecayMode = false;
			setReturnValue( "leptonicDecayMode" , askedLeptonicDecayMode );

			bosonDecayMode = findDecayMode( firstElectron , otherDecayProducts , leptonicDecayProducts );
			if ( bosonDecayMode <= 0 )
			{
				for ( unsigned int i_d = 0 ; i_d < firstElectron->getDaughters().size() ; ++i_d )
				{
					const EVENT::MCParticle *daughter = firstElectron->getDaughters()[ i_d ];
					if ( ( abs( daughter->getPDG() ) == 23 || abs( daughter->getPDG() ) == 25 ) && daughter != leptonicZ )
					{
						bosonDecayMode = findDecayMode( daughter , otherDecayProducts , leptonicDecayProducts );
					}
				}
			}
			if ( m_nTrueQuarks == 2 )
			{
				MCParticle* quark = dynamic_cast<MCParticle*>( TrueQuarkCollection->getElementAt( 0 ) );
				int quarkPDG = abs( quark->getPDG() );
				if ( quarkPDG == 1 ) bosonDecayMode = 6;
				if ( quarkPDG == 2 ) bosonDecayMode = 5;
				if ( quarkPDG == 3 ) bosonDecayMode = 4;
				if ( quarkPDG == 4 ) bosonDecayMode = 3;
				if ( quarkPDG == 5 ) bosonDecayMode = 2;
				if ( quarkPDG == 21 ) bosonDecayMode = 7;
			}
			if ( bosonDecayMode > 0 )
			{
				if ( m_includ_other && bosonDecayMode == 0 ) askedDecayMode = true;
				else if ( m_includ_WW && bosonDecayMode == 1 ) askedDecayMode = true;
				else if ( m_includ_bb && bosonDecayMode == 2 ) askedDecayMode = true;
				else if ( m_includ_cc && bosonDecayMode == 3 ) askedDecayMode = true;
				else if ( m_includ_ss && bosonDecayMode == 4 ) askedDecayMode = true;
				else if ( m_includ_uu && bosonDecayMode == 5 ) askedDecayMode = true;
				else if ( m_includ_dd && bosonDecayMode == 6 ) askedDecayMode = true;
				else if ( m_includ_gg && bosonDecayMode == 7 ) askedDecayMode = true;
				else if ( m_includ_ee && bosonDecayMode == 8 ) askedDecayMode = true;
				else if ( m_includ_mumu && bosonDecayMode == 9 ) askedDecayMode = true;
				else if ( m_includ_tautau && bosonDecayMode == 10 ) askedDecayMode = true;
				else if ( m_includ_nu1nu1 && bosonDecayMode == 11 ) askedDecayMode = true;
				else if ( m_includ_nu2nu2 && bosonDecayMode == 12 ) askedDecayMode = true;
				else if ( m_includ_nu3nu3 && bosonDecayMode == 13 ) askedDecayMode = true;
				else if ( m_includ_gammagamma && bosonDecayMode == 14 ) askedDecayMode = true;
				else  askedDecayMode = false;
				//setReturnValue( "ZHDecayMode" , askedDecayMode );
			}
			else
			{
				askedDecayMode = false;
			}
			setReturnValue( "ZHDecayMode" , askedDecayMode );
		}
		if ( m_fillRootTree )
		{
			h_ZHDecayMode->Fill( bosonDecayMode );
			n_ZHDecays += 1.0;
			h_ZLeptonicDecayMode->Fill( leptonicDecayMode );
			n_ZDecays += 1.0;
		}
		streamlog_out(MESSAGE) << "	Leptonic Decay Mode:		" << leptonicDecayMode << " (Di-Lepton True Invariant Mass = " << m_diLeptonTrueInvMass << " GeV , Reco Invariant Mass = " << m_diLeptonRecoInvMass << " GeV)" << std::endl;
		streamlog_out(MESSAGE) << "	Other (Boson) Decay Mode:	" << bosonDecayMode << " (Di-Quark Invariant Mass = " << m_diQuarkInvMass << " GeV)" << std::endl;
		streamlog_out(MESSAGE) << "	Event selected by number of IsoLeptons:		" << ( trueNIsoLeps ? "TRUE" : "FALSE" ) << std::endl;
		streamlog_out(MESSAGE) << "	Event selected by number of Jets:		" << ( trueNJets ? "TRUE" : "FALSE" ) << std::endl;
		streamlog_out(MESSAGE) << "	Event selected by Leptonic Decay Mode:		" << ( askedLeptonicDecayMode ? "TRUE" : "FALSE" ) << std::endl;
		streamlog_out(MESSAGE) << "	Event selected by Other (Boson) Decay Mode:	" << ( askedDecayMode ? "TRUE" : "FALSE" ) << std::endl;
		m_leptonicDecayMode = leptonicDecayMode;
		m_bosonDecayMode = bosonDecayMode;
	}
	catch( DataNotAvailableException &e )
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }
	if ( m_fillRootTree ) m_pTTree->Fill();

}

int ZHEvents::findDecayLeptonicMode( const EVENT::MCParticle *motherParticle , mcpVector &leptonicDecayProducts )
{
	int decayMode = -1;
	for ( unsigned int i_d1 = 0 ; i_d1 < motherParticle->getDaughters().size() ; ++i_d1 )
	{
		EVENT::MCParticle *firstDaughter = motherParticle->getDaughters()[ i_d1 ];
		for ( unsigned int i_type = 0 ; i_type < ( sizeof( leptonicPDGCodes ) / sizeof( *leptonicPDGCodes ) ) ; ++ i_type )
		{
			if ( firstDaughter->getPDG() == leptonicPDGCodes[ i_type ] )
			{
				for ( unsigned int i_d2 = 0 ; i_d2 < motherParticle->getDaughters().size() ; ++i_d2 )
				{
					if ( i_d1 != i_d2 )
					{
						EVENT::MCParticle *secondDaughter = motherParticle->getDaughters()[ i_d2 ];
						if ( firstDaughter->getPDG() == -1 * secondDaughter->getPDG() )
						{
							leptonicDecayProducts.push_back( firstDaughter );
							leptonicDecayProducts.push_back( secondDaughter );
							decayMode = i_type + 1 ;
						}
					}
				}
			}
		}
	}
	return decayMode;
}

int ZHEvents::findDecayMode( const EVENT::MCParticle *motherParticle , mcpVector &otherDecayProducts , mcpVector leptonicDecayProducts )
{
	int decayMode = 0;
	for ( unsigned int i_d1 = 0 ; i_d1 < motherParticle->getDaughters().size() ; ++i_d1 )
	{
		EVENT::MCParticle *firstDaughter = motherParticle->getDaughters()[ i_d1 ];
		bool isFirstDaughterInLeptonicMode = false;
		for ( unsigned int i_lep = 0 ; i_lep < leptonicDecayProducts.size() ; ++i_lep )
		{
			if ( firstDaughter == leptonicDecayProducts[ i_lep ] ) isFirstDaughterInLeptonicMode = true;
		}
		if ( !isFirstDaughterInLeptonicMode )
		{
			for ( unsigned int i_type = 0 ; i_type < ( sizeof( PDGCodes ) / sizeof( *PDGCodes ) ) ; ++ i_type )
			{
				if ( firstDaughter->getPDG() == PDGCodes[ i_type ] )
				{
					for ( unsigned int i_d2 = 0 ; i_d2 < motherParticle->getDaughters().size() ; ++i_d2 )
					{
						if ( i_d1 != i_d2 )
						{
							EVENT::MCParticle *secondDaughter = motherParticle->getDaughters()[ i_d2 ];
							bool isSecondDaughterInLeptonicMode = false;
							for ( unsigned int i_lep = 0 ; i_lep < leptonicDecayProducts.size() ; ++i_lep )
							{
								if ( secondDaughter == leptonicDecayProducts[ i_lep ] ) isSecondDaughterInLeptonicMode = true;
							}
							if ( !isSecondDaughterInLeptonicMode )
							{
								if ( abs( firstDaughter->getPDG() ) == 24 && firstDaughter->getPDG() == -1 * secondDaughter->getPDG() )
								{
									decayMode = i_type + 1 ;
									otherDecayProducts.push_back( firstDaughter );
									otherDecayProducts.push_back( secondDaughter );
								}
								else if ( firstDaughter->getPDG() <= 16 && firstDaughter->getPDG() == -1 * secondDaughter->getPDG() )
								{
									decayMode = i_type + 1 ;
									otherDecayProducts.push_back( firstDaughter );
									otherDecayProducts.push_back( secondDaughter );
								}
								else if ( firstDaughter->getPDG() > 16 && firstDaughter->getPDG() == secondDaughter->getPDG() )
								{
									decayMode = i_type + 1 ;
									otherDecayProducts.push_back( firstDaughter );
									otherDecayProducts.push_back( secondDaughter );
								}
							}
						}
					}
				}
			}
		}
	}
	return decayMode;
}

void ZHEvents::check()
{
}

void ZHEvents::end()
{
	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		//h_ZHDecayMode->Scale( 100.0 / n_ZHDecays );
		h_ZHDecayMode->Write();
		//h_ZLeptonicDecayMode->Scale( 100.0 / n_ZDecays );
		h_ZLeptonicDecayMode->Write();
		m_pTFile->Close();
		delete m_pTFile;
	}
}
