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
// -----------------------------------------
// This is a wrapper for iS3D with the JETSCAPE framework
// -----------------------------------------

#include "JetScapeLogger.h"
#include "iS3DWrapper.h"

#include <string>

using namespace Jetscape;

iS3DWrapper::iS3DWrapper() {
    SetId("iS3D");
    //iS3D_ptr_ = nullptr;
}

iS3DWrapper::~iS3DWrapper() {
}

void iS3DWrapper::InitTask() {
    INFO << "Initialize a particle sampler (iS3D)";

    /*
    iS3D_xml_ = xml_->FirstChildElement("iS3D");
    if (!iS3D_xml_) {
        WARN << "No XML section for iS3D! Please check the input file~";
        exit(-1);
    }
    string input_file = (iS3D_xml_->FirstChildElement("iS3D_input_file")->GetText());
    string working_path = (iS3D_xml_->FirstChildElement("iS3D_working_path")->GetText());
    int hydro_mode;
    iS3D_xml_->FirstChildElement("hydro_mode")->QueryIntText(&hydro_mode);

    int number_of_repeated_sampling;
    iS3D_xml_->FirstChildElement("number_of_repeated_sampling")->QueryIntText(&number_of_repeated_sampling);

    int flag_perform_decays;
    iS3D_xml_->FirstChildElement("Perform_resonance_decays")->QueryIntText(&flag_perform_decays);

    iS3D_ptr_ = new iS3D(working_path);
    iS3D_ptr_->paraRdr_ptr->readFromFile(input_file);

    // overwrite some parameters
    iS3D_ptr_->paraRdr_ptr->setVal("hydro_mode", hydro_mode);
    iS3D_ptr_->paraRdr_ptr->setVal("output_samples_into_files", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("use_OSCAR_format", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("use_gzip_format", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("store_samples_in_memory", 1);
    iS3D_ptr_->paraRdr_ptr->setVal("number_of_repeated_sampling", number_of_repeated_sampling);
    iS3D_ptr_->paraRdr_ptr->setVal("perform_decays", flag_perform_decays);

    // set default parameters
    iS3D_ptr_->paraRdr_ptr->setVal("turn_on_shear", 1);
    iS3D_ptr_->paraRdr_ptr->setVal("turn_on_bulk", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("turn_on_rhob", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("turn_on_diff", 0);

    iS3D_ptr_->paraRdr_ptr->setVal("include_deltaf_shear", 1);
    iS3D_ptr_->paraRdr_ptr->setVal("include_deltaf_bulk", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("bulk_deltaf_kind", 1);
    iS3D_ptr_->paraRdr_ptr->setVal("include_deltaf_diffusion", 0);

    iS3D_ptr_->paraRdr_ptr->setVal("restrict_deltaf", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("deltaf_max_ratio", 1.0);
    iS3D_ptr_->paraRdr_ptr->setVal("f0_is_not_small", 1);

    iS3D_ptr_->paraRdr_ptr->setVal("calculate_vn", 0);
    iS3D_ptr_->paraRdr_ptr->setVal("MC_sampling", 2);

    iS3D_ptr_->paraRdr_ptr->setVal("sample_upto_desired_particle_number", 0);
    iS3D_ptr_->paraRdr_ptr->echo();

    */
}

void iS3DWrapper::Exec() {

  /*
    int status = iS3D_ptr_->read_in_FO_surface();
    if (status != 0) {
        WARN << "Some errors happened in reading in the hyper-surface";
        exit(-1);
    }

    auto random_seed = (*GetMt19937Generator())();  // get random seed
    iS3D_ptr_->set_random_seed(random_seed);
    VERBOSE(2) << "Random seed used for the iS3D module" << random_seed;

    status = iS3D_ptr_->generate_samples();
    if (status != 0) {
        WARN << "Some errors happened in generating particle samples";
        exit(-1);
    }
    PassHadronListToJetscape();

    */

  // argument 0 : read surface from memory
  // argument 1 : read fo surface from file
  iS3D.run_particlization(1);
  PassHadronListToJetscape();
}

void iS3DWrapper::Clear() {
    VERBOSE(2) << "Finish iS3D particle sampling";
    //if (iS3D_ptr_ != nullptr) delete iS3D_ptr_;
}

void iS3DWrapper::PassHadronListToJetscape()
{
  /*
  unsigned int nev = iS3D_ptr_->get_number_of_sampled_events();
  VERBOSE(2) << "Passing all sampled hadrons to the JETSCAPE framework";
  VERBOSE(4) << "number of events to pass : " << nev;
  for (unsigned int iev = 0; iev < nev; iev++) {
  std::vector<shared_ptr<Hadron>> hadrons;
  unsigned int nparticles = (iS3D_ptr_->get_number_of_particles(iev));
  VERBOSE(4) << "event " << iev << ": number of particles = "<< nparticles;
  for (unsigned int ipart = 0; ipart < nparticles; ipart++) {
  iSS_Hadron current_hadron = (iS3D_ptr_->get_hadron(iev, ipart));
  int hadron_label = 0;
  int hadron_status = -1;
  int hadron_id = current_hadron.pid;
  //int hadron_id = 1;   // just for testing need to be changed to the line above
  double hadron_mass = current_hadron.mass;
  FourVector hadron_p(current_hadron.px, current_hadron.py, current_hadron.pz, current_hadron.E);
  FourVector hadron_x(current_hadron.x, current_hadron.y, current_hadron.z, current_hadron.t);

  // create a JETSCAPE Hadron
  hadrons.push_back(make_shared<Hadron>(hadron_label,hadron_id,hadron_status,hadron_p,hadron_x, hadron_mass));
  //Hadron* jetscape_hadron = new Hadron(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass);
  //(*Hadron_list_)[iev]->push_back(*jetscape_hadron);
}
Hadron_list_.push_back(hadrons);
}
VERBOSE(4) << "JETSCAPE received " << Hadron_list_.size() << " events.";
for (unsigned int iev = 0; iev < Hadron_list_.size(); iev++) {
VERBOSE(4) << "In event " << iev << " JETSCAPE received " << Hadron_list_.at(iev).size() << " particles.";
}
*/


  VERBOSE(2) << "Passing all sampled hadrons to the JETSCAPE framework";
  std::vector<shared_ptr<Hadron>> hadrons;
  unsigned int nparticles = ( iS3D.final_particles_.size() );
  VERBOSE(4) << "number of particles = " << nparticles;
  for (unsigned int ipart = 0; ipart < nparticles; ipart++)
  {
    Sampled_Particle current_hadron = iS3D.final_particles_[ipart];
    int hadron_label = 0;
    int hadron_status = -1;
    int hadron_id = current_hadron.mcID;
    double hadron_mass = current_hadron.mass;
    FourVector hadron_p(current_hadron.px, current_hadron.py, current_hadron.pz, current_hadron.E);
    FourVector hadron_x(current_hadron.x, current_hadron.y, current_hadron.z, current_hadron.t);

    // create a JETSCAPE Hadron
    hadrons.push_back(make_shared<Hadron>(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass));
    //Hadron* jetscape_hadron = new Hadron(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass);
    //(*Hadron_list_)[iev]->push_back(*jetscape_hadron);
  }
  Hadron_list_.push_back(hadrons);

  VERBOSE(4) << "JETSCAPE received " << Hadron_list_.size() << " events.";
  for (unsigned int iev = 0; iev < Hadron_list_.size(); iev++)
  {
    VERBOSE(4) << "In event " << iev << " JETSCAPE received " << Hadron_list_.at(iev).size() << " particles.";
  }
}

void iS3DWrapper::WriteTask(weak_ptr<JetScapeWriter> w)
{
  /*
  VERBOSE(4)<<"In iS3DWrapper::WriteTask";
  auto f = w.lock();
  if ( !f ) return;

  f->WriteComment("JetScape module: "+GetId());
  if(Hadron_list_.size()>0) {
    f->WriteComment("Final State Bulk Hadrons");
    for(unsigned int j=0; j<Hadron_list_.size(); j++){
      vector<shared_ptr<Hadron>> hadVec = Hadron_list_.at(j);
      for(unsigned int i=0; i<hadVec.size(); i++) {
	f->WriteWhiteSpace("["+to_string(i)+"] H");
	f->Write(hadVec.at(i));
      }
    }
  } else {
    f->WriteComment("There are no bulk Hadrons");
  }
 */
}