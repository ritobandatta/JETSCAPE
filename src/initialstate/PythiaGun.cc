/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

// Create a pythia collision at a specified point and return the two inital hard partons

#include "PythiaGun.h"
#include <sstream>
#include <iostream>
#include <fstream>
#define MAGENTA "\033[35m"

using namespace std;

int pthat_bins[33][2] = {{1500,2510}, {800,1500}, {450, 800}, {350, 450}, {300, 350}, {280, 300}, {260, 280}, {240, 260}, {220, 240}, {210, 220}, {200, 210}, {190, 200}, {180,190}, {170, 180}, {160, 170}, {150, 160}, {140, 150}, {130, 140}, {120, 130}, {110, 120},  {100, 110}, {90,100},  {80,90}, {70,80}, {60, 70}, {55, 60}, {50, 55},  {45, 50}, {40, 45}, {35, 40}, {30, 35}, {25, 30}, {20,25}};

// Register the module with the base class
RegisterJetScapeModule<PythiaGun> PythiaGun::reg("PythiaGun");

PythiaGun::~PythiaGun() { VERBOSE(8); }

void PythiaGun::InitTask() {

  JSDEBUG << "Initialize PythiaGun";
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = off");
  readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetInfo()) {
    readString("Init:showProcesses = on");
    readString("Init:showChangedSettings = on");
    readString("Init:showMultipartonInteractions = on");
    readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  std::string s = GetXMLElementText({"Hard", "PythiaGun", "name"});
  SetId(s);
  // cout << s << endl;

  //other Pythia settings
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = off");
  readString("PartonLevel:ISR = on");
  readString("PartonLevel:MPI = on");
  //readString("PartonLevel:FSR = on");
  readString("PromptPhoton:all=on");
  readString("WeakSingleBoson:all=off");
  readString("WeakDoubleBoson:all=off");

  pTHatMin = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMin"});
  pTHatMax = GetXMLElementDouble({"Hard", "PythiaGun", "pTHatMax"});

  if(pTHatMin < 0.01){ //assuming low bin where softQCD should be used
    //running softQCD - inelastic nondiffrative (min-bias)
    readString("HardQCD:all = off");
    readString("SoftQCD:nonDiffractive = on");
    softQCD = true;
  }
  else{ //running normal hardQCD
    readString("HardQCD:all = on"); // will repeat this line in the xml for demonstration
    //  readString("HardQCD:gg2ccbar = on"); // switch on heavy quark channel
    //readString("HardQCD:qqbar2ccbar = on");
    numbf.str("PhaseSpace:pTHatMin = "); numbf << pTHatMin; readString(numbf.str());
    numbf.str("PhaseSpace:pTHatMax = "); numbf << pTHatMax; readString(numbf.str());
	  softQCD = false;
  }

  // SC: read flag for FSR
  FSR_on = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
  if (FSR_on)
    readString("PartonLevel:FSR = on");
  else
    readString("PartonLevel:FSR = off");

  JSINFO << MAGENTA << "Pythia Gun with FSR_on: " << FSR_on;
  JSINFO << MAGENTA << "Pythia Gun with " << pTHatMin << " < pTHat < " << pTHatMax;

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if (RandomXmlDescription) {
    tinyxml2::XMLElement *xmle =
        RandomXmlDescription->FirstChildElement("seed");
    if (!xmle)
      throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    JSWARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) << "Seeding pythia to " << seed;
  numbi << seed;
  readString(numbi.str());

  // Species
  readString("Beams:idA = 2212");
  readString("Beams:idB = 2212");

  // Energy
  eCM = GetXMLElementDouble({"Hard", "PythiaGun", "eCM"});
  numbf.str("Beams:eCM = ");
  numbf << eCM;
  readString(numbf.str());

  // Reading vir_factor from xml for MATTER
  vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});
  initial_virtuality_pT = GetXMLElementInt({"Eloss", "Matter", "initial_virtuality_pT"});
  if(vir_factor < rounding_error) {
    JSWARN << "vir_factor should not be zero or negative";
    exit(1);
  }

  softMomentumCutoff = GetXMLElementDouble({"Hard", "PythiaGun", "softMomentumCutoff"});

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "PythiaGun", "LinesToRead"}, false);
  int i = 0;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  // And initialize
  if (!init()) { // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }

    std::ofstream sigma_printer;
    sigma_printer.open(printer, std::ios::trunc);


}

void PythiaGun::Exec() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  pTHatMin=pthat_bins[GetCurrentEvent()/100][0];
  pTHatMax=pthat_bins[GetCurrentEvent()/100][1];
	
  bool flag62 = false;
  vector<Pythia8::Particle> p62;

  // sort by pt
  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };

  do {
    next();
    p62.clear();
      if (!printer.empty()){
            std::ofstream sigma_printer;
            sigma_printer.open(printer, std::ios::out | std::ios::app);

            sigma_printer << "sigma = " << GetSigmaGen() << " Err =  " << GetSigmaErr() << endl ;
            //sigma_printer.close();

//      JSINFO << BOLDYELLOW << " sigma = " << GetSigmaGen() << " sigma err = " << GetSigmaErr() << " printer = " << printer << " is " << sigma_printer.is_open() ;
    };

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for (int parid = 0; parid < event.size(); parid++) {
      if (parid < 3)
        continue; // 0, 1, 2: total event and beams
      Pythia8::Particle &particle = event[parid];

      //replacing diquarks with antiquarks (and anti-dq's with quarks)
      //the id is set to the heaviest quark in the diquark (except down quark)
      //this technically violates baryon number conservation over the entire event
      //also can violate electric charge conservation
      if( (std::abs(particle.id()) > 1100) && (std::abs(particle.id()) < 6000) && ((std::abs(particle.id())/10)%10 == 0) ){
        if(particle.id() > 0){particle.id( -1*particle.id()/1000 );}
        else{particle.id( particle.id()/1000 );}
      }

      if (!FSR_on) {
        // only accept particles after MPI
        if (particle.status() != 62)
          continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;

        // reject rare cases of very soft particles that don't have enough e to get
        // reasonable virtuality
        if (initial_virtuality_pT && (particle.pT() < softMomentumCutoff)) {
          // this cutoff was 1.0/sqrt(vir_factor) in versions < 3.6
          continue;
        } else if(!initial_virtuality_pT && (particle.pAbs() < softMomentumCutoff)) {
          continue;
        }

        //if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
        // accept
      } else { // FSR_on true: use Pythia vacuum shower instead of MATTER
        if (!particle.isFinal())
          continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;
      }
      p62.push_back(particle);
    }

    // if you want at least 2
    if (p62.size() < 2)
      continue;
    //if ( p62.size() < 1 ) continue;

    // Now have all candidates, sort them
    // sort by pt
    std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    //skipping event if softQCD is on & pThat exceeds max (where next bin is HardQCD with this as pThatmin)
    if(softQCD && (info.pTHat() >= pTHatMax)){continue;}

    flag62 = true;

  } while (!flag62);

  double p[4], xLoc[4];

  // This location should come from an initial state
  for (int i = 0; i <= 3; i++) {
    xLoc[i] = 0.0;
  };

  // // Roll for a starting point
  // // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

  if (!ini) {
    VERBOSE(1) << "No initial state module, setting the starting location to "
                  "0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    double x, y;
    ini->SampleABinaryCollisionPoint(x, y);
    xLoc[1] = x;
    xLoc[2] = y;
  }

  // Loop through particles

  // Only top two
  //for(int np = 0; np<2; ++np){

  // Accept them all

  int hCounter = 0;
  for (int np = 0; np < p62.size(); ++np) {
    Pythia8::Particle &particle = p62.at(np);

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << ", pT = " << particle.pT() << ", eta = " << particle.eta()
               << ", phi = " << particle.phi() << ", e = " << particle.e();

    VERBOSE(7) << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    auto ptn = make_shared<Parton>(0, particle.id(), 0, particle.pT(), particle.eta(),particle.phi(), particle.e(), xLoc);
    ptn->set_color(particle.col());
    ptn->set_anti_color(particle.acol());
    ptn->set_max_color(1000 * (np + 1));
    AddParton(ptn);
  }

  VERBOSE(8) << GetNHardPartons();
}
