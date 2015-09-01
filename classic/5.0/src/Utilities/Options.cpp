
namespace Options {
    // The global program options.
    bool echo = false;
    bool info = true;
    bool csrDump = false;
    bool ppdebug = false;

    // If true create symmetric distribution
    bool cZero = false;

    bool enableHDF5 = true;
    bool asciidump = false;

    // If the distance of a particle to bunch mass larger than remotePartDel times of the rms size of the bunch in any dimension,
    // the particle will be deleted artifically to hold the accuracy of space charge calculation. The default setting of -1 stands for no deletion.
    int remotePartDel = -1;

  double beamHaloBoundary = 0;

}
