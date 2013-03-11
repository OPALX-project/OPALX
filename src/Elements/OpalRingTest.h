
class OpalRingTest {
public:
  RunTests();

  OpalRing* setupRingHardCoded() {
      using Physics::pi;
      double phi_init = -pi;
      double phi_total = phi_init;
      lattice_rInit_m = 2500.;
      for (size_t i = 0; i < 8; ++i) {
          Vector_t position(0.0);
          position(0) = sin(phi_total)*lattice_rInit_m;;
          position(2) = cos(phi_total)*lattice_rInit_m;;
          Vector_t orientation(0.0);
          orientation(0) = cos(phi_total);
          orientation(2) = -sin(phi_total);
          Component* bend = new SBend3D("Name");
          OpalRingSection* section = new OpalRingSection();
          section->setComponent(bend);
          section->setStartPosition(position);
          section->setStartNormal(orientation);
          section_list_m.push_back(section); // compvec, start_phi, end_phi - what does start/end mean?
          if (i > 0) {
              section_list_m[i-1]->setEndPosition(position);
              section_list_m[i-1]->setEndNormal(orientation);
          }
          std::cerr << "SETUP " << i << " pos: " << position << " n: " << orientation << " phi: " << phi_total << " ITEM: " << section_list_m[i] << std::endl;
          phi_total += pi/4.;
      }
      if (fabs(phi_total-phi_init-2*pi) > 1e-3) {
          throw OpalException("OpalRing::setupSectionHardCoded()",
                              "Ring geometry does not sum to 2*PI");
      }
  }

  TestOpalRing() {
    section_list_m.back()->setEndPosition(section_list_m[0]->getStartPosition());
    section_list_m.back()->setEndNormal(section_list_m[0]->getStartNormal());
    current_section_m = section_list_m.begin();
    // expect either (phi_total-phi_init) is 2.*pi within tolerance
    for (size_t i = 0; i < 8; ++i) {
        double test_phi = pi/8.+i*pi/4.;
        Vector_t position(0.);
        position(0) = sin(test_phi)*lattice_rInit_m;
        position(2) = cos(test_phi)*lattice_rInit_m;
        OpalRingSection *section = getSectionAt(position);
        std::cerr << "TEST " << i << " pos: " << position << " item: " << section << std::endl;
    }
}

};

