#include "MultiBunchDump.h"

#include <boost/filesystem.hpp>

#include <iomanip>

#include "AbstractObjects/OpalData.h"
#include "Utilities/Timer.h"

MultiBunchDump::MultiBunchDump()
    : nBins_m(-1)
    , fbase_m(OpalData::getInstance()->getInputBasename())
    , fext_m(".smb")
{ }


void MultiBunchDump::writeHeader(const std::string& fname) const {
    
    if ( boost::filesystem::exists(fname) ) {
        return;
    }
    
    std::ofstream out;
    this->open_m(out, fname);
    
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    std::string indent("        ");
    
    out << "SDDS1" << std::endl;
    out << "&description\n"
        << indent << "text=\"Multi Bunch Statistics data '" << OpalData::getInstance()->getInputFn()
        << "' " << dateStr << " " << timeStr << "\",\n"
        << indent << "contents=\"multi bunch stat parameters\"\n"
        << "&end\n";
    out << "&parameter\n"
        << indent << "name=processors,\n"
        << indent << "type=long,\n"
        << indent << "description=\"Number of Cores used\"\n"
        << "&end\n";
    out << "&parameter\n"
        << indent << "name=revision,\n"
        << indent << "type=string,\n"
        << indent << "description=\"git revision of opal\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=t,\n"
        << indent << "type=double,\n"
        << indent << "units=ns,\n"
        << indent << "description=\"1 Time\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=numParticles,\n"
        << indent << "type=long,\n"
        << indent << "units=1,\n"
        << indent << "description=\"2 Number of Macro Particles\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=energy,\n"
        << indent << "type=double,\n"
        << indent << "units=MeV,\n"
        << indent << "description=\"3 Mean Bunch Energy\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=dE,\n"
        << indent << "type=double,\n"
        << indent << "units=MeV,\n"
        << indent << "description=\"16 energy spread of the beam\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_x,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"4 RMS Beamsize in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_y,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"5 RMS Beamsize in y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=rms_s,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"6 RMS Beamsize in s\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=emit_x,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"7 Normalized Emittance x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=emit_y,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"8 Normalized Emittance y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=emit_s,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"9 Normalized Emittance s\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=mean_x,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"10 Mean Beam Position in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=mean_y,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"11 Mean Beam Position in y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=mean_s,\n"
        << indent << "type=double,\n"
        << indent << "units=m,\n"
        << indent << "description=\"12 Mean Beam Position in s\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=halo_x,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"13 Halo in x\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=halo_y,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"14 Halo in y\"\n"
        << "&end\n";
    out << "&column\n"
        << indent << "name=halo_z,\n"
        << indent << "type=double,\n"
        << indent << "units=1,\n"
        << indent << "description=\"15 Halo in z\"\n"
        << "&end\n";
    
    close_m(out);
}


void MultiBunchDump::writeData(const beaminfo_t& binfo, int bin) {
    
    std::string fname = fbase_m + "-bin-" + std::to_string(bin) + fext_m;
    
    writeHeader(fname);
    
    std::ofstream out;
    open_m(out, fname);
    
    unsigned int pwi = 10;
    
    out << binfo.time       << std::setw(pwi) << "\t"
        << binfo.nParticles << std::setw(pwi) << "\t"
        << binfo.ekin       << std::setw(pwi) << "\t"
        << binfo.dEkin      << std::setw(pwi) << "\t"
        << binfo.rms[0]     << std::setw(pwi) << "\t"
        << binfo.rms[1]     << std::setw(pwi) << "\t"
        << binfo.rms[2]     << std::setw(pwi) << "\t"
        << binfo.emit[0]    << std::setw(pwi) << "\t"
        << binfo.emit[1]    << std::setw(pwi) << "\t"
        << binfo.emit[2]    << std::setw(pwi) << "\t"
        << binfo.mean[0]    << std::setw(pwi) << "\t"
        << binfo.mean[1]    << std::setw(pwi) << "\t"
        << binfo.mean[2]    << std::setw(pwi) << "\t"
        << binfo.halo[0]    << std::setw(pwi) << "\t"
        << binfo.halo[1]    << std::setw(pwi) << "\t"
        << binfo.halo[2]    << std::endl;
    
    close_m(out);
}


void MultiBunchDump::open_m(std::ofstream& out,
                            const std::string& fname) const
{
    std::ios::openmode mode = std::ios::out;
    if ( boost::filesystem::exists(fname) ) {
        mode = std::ios::app;
    }
    
    out.open(fname.c_str(), mode);
    out.precision(15);
    out.setf(std::ios::scientific, std::ios::floatfield);
}


void MultiBunchDump::close_m(std::ofstream& out) const {
    if ( out.is_open() ) {
        out.close();
    }
}
