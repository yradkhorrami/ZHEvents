#ifndef ZHEvents_h_1
#define ZHEvents_h_1

#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "lcio.h"
#include "EVENT/LCStrVec.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "marlin/VerbosityLevels.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "TFile.h"
#include "TH1I.h"
#include "TTree.h"
#include "TLorentzVector.h"

class TFile;
class TH1F;
class TH1I;
class TTree;

using namespace lcio ;
using namespace marlin ;

class ZHEvents : public Processor
{
public:
	virtual Processor *newProcessor()
	{
		return new ZHEvents;
	}
	ZHEvents();
	virtual ~ZHEvents() = default;
	ZHEvents( const ZHEvents& ) = delete;
	ZHEvents &operator = ( const ZHEvents& ) = delete;
	typedef std::vector<EVENT::MCParticle*> mcpVector;
	virtual void init();
	virtual void Clear();
	virtual void processRunHeader();
	virtual void processEvent( EVENT::LCEvent *pLCEvent );
	int findDecayLeptonicMode( const EVENT::MCParticle *motherParticle , mcpVector &leptonicDecayProducts );
	int findDecayMode( const EVENT::MCParticle *motherParticle , mcpVector &otherDecayProducts , mcpVector leptonicDecayProducts );
 	virtual void check();
	virtual void end();
private:

	typedef std::vector<int>		IntVector;
	typedef std::vector<double>		DoubleVector;
	typedef std::vector<float>		FloatVector;

	std::string				m_mcParticleCollection{};
	std::string				m_inputJetCollection{};
	std::string				m_inputIsoLepCollection{};
	std::string				m_inputTrueQuarkCollection{};
	std::string				m_inputTrueIsoLepCollection{};
	std::string				m_rootFile{};

	int					m_nRun;
	int					m_nEvt;
	int					m_nRunSum;
	int					m_nEvtSum;
	bool					m_cheatDecayMode = true;
	bool					m_fillRootTree = true;

	bool					m_includ_bb = true;
	bool					m_includ_cc = true;
	bool					m_includ_ss = true;
	bool					m_includ_uu = true;
	bool					m_includ_dd = true;
	bool					m_includ_gg = true;
	bool					m_includ_ee = true;
	bool					m_includ_mumu = true;
	bool					m_includ_tautau = true;
	bool					m_includ_nu1nu1 = true;
	bool					m_includ_nu2nu2 = true;
	bool					m_includ_nu3nu3 = true;
	bool					m_includ_gammagamma = true;
	bool					m_includ_WW = true;
	bool					m_includ_ZZ = true;
	bool					m_includ_HH = true;
	bool					m_includ_other = true;
	bool					m_includZe1e1 = true;
	bool					m_includZe2e2 = true;
	bool					m_includZe3e3 = true;
	int					m_nJets = 0;
	int					m_nIsoLeps = 0;
	int					m_nRecoIsoLeptons;
	int					m_nRecoJets;
	int					m_nTrueIsoLeptons;
	int					m_nTrueQuarks;
	int					m_leptonicDecayMode;
	int					m_bosonDecayMode;
	float					m_diLeptonInvMass;
	float					m_diQuarkInvMass;
	float					m_isoLepInvMassCutMin;
	float					m_isoLepInvMassCutMax;

	double					PDGCodes[ 14 ]{ 24 , 5 , 4 , 3 , 2 , 1 , 21 , 11 , 13 , 15 , 12 , 14 , 16 , 22 };
	double					leptonicPDGCodes[ 3 ]{ 11 , 13 , 15 };
	TFile					*m_pTFile{};
	TTree					*m_pTTree{};
	TH1F					*h_ZHDecayMode{};
	float					n_ZHDecays = 0.0;
	TH1F					*h_ZLeptonicDecayMode{};
	float					n_ZDecays = 0.0;
};
#endif
